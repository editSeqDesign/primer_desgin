import pandas as pd
from Bio import SeqIO 
import os,sys
import sgRNA_utils.sgRNA_primer_util as su
import module.sequencing_primer as sr
import module.from_gbFile_to_seq as fq
import module.parser_input_to_df as pf
import module.product_and_decorate_editingSeq as p_d_seq
import module.order as order
import warnings   
warnings.filterwarnings('ignore')  
import configparser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from sgRNA_utils.sgRNA_primer_config import config 
from loguru import logger

#uha_dha_primer
def extract_uha_dha_primer(info_input_df, sgRNA):

    primer_template_for_u_d_ha_df = p_d_seq.create_primer_template(info_input_df,sgRNA)
    uha_dha_primer_df = p_d_seq.design_primer(primer_template_for_u_d_ha_df,'Name_Region','primer_template','u_d')
    uha_dha_primer_df = uha_dha_primer_df.rename(columns={'Region':'id'})
    uha_dha_primer_df = uha_dha_primer_df.join(uha_dha_primer_df.id.str.split(';',expand=True).rename(columns = {0:'Name',1:'Region',2:'Type'})).drop(columns='id')
    
    return uha_dha_primer_df

#uha_dha
def extract_uha_dha(info_input_df,uha_dha_primer_df,sgRNA):
    #整合进突变信息
    info_df = info_input_df[['name','region','seq_altered','type','ref','strand']].rename(columns={'name':'Name','region':'Region'})
    info_df.seq_altered.fillna('',inplace=True)
    uha_dha_info_primer_df = pd.merge(info_df,uha_dha_primer_df,on=['Name','Region'])

    #提取源生同源臂
    uha_dha_df = p_d_seq.create_uha_dha_df(uha_dha_primer_df) 
    #合并突变信息
    uha_dha_df = pd.merge(uha_dha_df,info_df,on=['Name','Region'])
    uha_dha_sgRNA_df = pd.merge(uha_dha_df,sgRNA,on=['Name','Region'],how='inner')

    return uha_dha_info_primer_df, uha_dha_df, uha_dha_sgRNA_df,info_df

#one plasmid pcr primer  
def one_plasmid_system_pcr_design_primer(gb_path,info_df,uha_dha_sgRNA_df,uha_dha_info_primer_df,uha_dha_primer_df,enzyme_df,enzyme_name):
     #创建新的质粒
    uha_dha_sgRNA_df, promoter_terminator_up_promoter_seq, promoter_terminator_down_terminator_seq, type_kind  = p_d_seq.create_new_plasmid(gb_path, uha_dha_sgRNA_df.copy(), ccdb_label='ccdB', promoter_terminator_label='gRNA', n_20_label='N20')

    #设计sgRNA、ccdb质粒引物
    n20up_primer_template = uha_dha_sgRNA_df[['Name','Region','n20_up_template','Target sequence','Rev Target sequence']]
    n20up_primer_template['Region'] = n20up_primer_template['Name'] +';'+ n20up_primer_template['Region']
    n20up_primer_df = p_d_seq.design_primer(n20up_primer_template,'Region','n20_up_template','sgRNA')
    n20up_primer_df = pd.merge(n20up_primer_template[['Region','Target sequence','Rev Target sequence']],n20up_primer_df,on=['Region'],how='inner')


    n20down_primer_template = uha_dha_sgRNA_df[['Name','Region','n20_down_template','Target sequence','Rev Target sequence']]
    n20down_primer_template['Region'] = n20down_primer_template['Name'] +';'+ n20down_primer_template['Region']
    n20down_primer_df = p_d_seq.design_primer(n20down_primer_template,'Region','n20_down_template','sgRNA')
    n20down_primer_df = pd.merge(n20down_primer_template[['Region','Target sequence','Rev Target sequence']],n20down_primer_df,on=['Region'],how='inner')

    seq_altered_primer_template = info_df[info_df.seq_altered.apply(lambda x:len(x)>120)][['Name','Region','seq_altered']]
    seq_altered_primer_template['Region'] = seq_altered_primer_template['Name'] + seq_altered_primer_template['Region']
    seq_altered_primer_df = p_d_seq.design_primer(seq_altered_primer_template,'Region','seq_altered',stype='seq_altered')

    #给同源臂引物加接头:uha取promoter_terminator_down_terminator_seq尾部反义4bp，dha取头promoter_terminator_up_promoter_seq正义4bp
    uha_dha_primer_df = p_d_seq.add_joint_sgRNA_primer(uha_dha_info_primer_df,enzyme_df,enzyme_name,promoter_terminator_down_terminator_seq, promoter_terminator_up_promoter_seq, stype = 'u_d_primer_joint')

    #ccdb、sgrna质粒引物加接头
    n20up_primer_df = p_d_seq.add_joint_sgRNA_primer(n20up_primer_df,enzyme_df,enzyme_name,'',stype='n20up_primer_joint')
    n20down_primer_df = p_d_seq.add_joint_sgRNA_primer(n20down_primer_df,enzyme_df,enzyme_name,'',stype='n20down_primer_joint')

    #seq_altered_primer加接头
    seq_altered_primer_df = p_d_seq.add_joint_sgRNA_primer(seq_altered_primer_df,enzyme_df,enzyme_name,'',stype='seq_altered_primer_joint')

    #分别提取，修饰后的uha、dha引物
    in_col = ['Name','Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']
    ou_col = ['Name','Region',"u_primer_f_seq_(5'-3')","u_primer_r_seq_(5'-3')",'UHA','UHA_size',"d_primer_f_seq_(5'-3')","d_primer_r_seq_(5'-3')",'DHA','DHA_size']
    uha_primer_df, dha_primer_df = p_d_seq.create_uha_dha_primer_df(uha_dha_primer_df, in_col, ou_col)

    n20down_primer_p_df = n20down_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]
    n20up_primer_p_df = n20up_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]

    return uha_dha_sgRNA_df,uha_primer_df, dha_primer_df, n20down_primer_p_df, n20up_primer_p_df, seq_altered_primer_df,type_kind

#one plasmid sequencing primer
def one_plasmid_system_sequencing_design_primer(type_kind,uha_dha_sgRNA_df):
    print(type_kind,"1:代表启动子上游的序列大于600bp,同时终止子下游大于600bp","2:代表启动子上游的序列小于600bp","3:代表终止子上游序列小于600bp")

    if type_kind == 3:
        #载体测序
        sequencing_primer_df=uha_dha_sgRNA_df[['Name', 'Region', 'UHA', 'UHA_size', 'DHA', 'DHA_size',
                                    'plasmid', 'promoter_N20_terminator',
                                    'promoter_N20_terminator_up', 'promoter_N20_terminator_down','seq_altered']]
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['UHA']+sequencing_primer_df['seq_altered']+sequencing_primer_df['DHA']+sequencing_primer_df['promoter_N20_terminator_up'] + sequencing_primer_df['promoter_N20_terminator']
        sequencing_primer_df['Region'] = sequencing_primer_df['Name']+';'+sequencing_primer_df['Region']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region')

    elif type_kind == 1:
        #载体测序
        #uha_dha测序
        sequencing_primer_df=uha_dha_sgRNA_df[['Name', 'Region', 'UHA', 'UHA_size', 'DHA', 'DHA_size',
                                    'plasmid', 'promoter_N20_terminator',
                                    'promoter_N20_terminator_up', 'promoter_N20_terminator_down','seq_altered']]
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['UHA'] + sequencing_primer_df['seq_altered'] + sequencing_primer_df['DHA']
        sequencing_primer_df['Region'] = sequencing_primer_df['Name']+';'+sequencing_primer_df['Region']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df1 = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region')
        #promoter_N20_terminator测序
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['promoter_N20_terminator']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df2 = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region')
        plasmid_sequencing_primer_df = p_d_seq.merge_sequencing_result(plasmid_sequencing_primer_df1, plasmid_sequencing_primer_df2)

    elif type_kind == 2:
        #载体测序
        #uha_dha测序
        sequencing_primer_df=uha_dha_sgRNA_df[[ 'Name', 'Region', 'UHA', 'UHA_size', 'DHA', 'DHA_size',
                                                'plasmid', 'promoter_N20_terminator',
                                                'promoter_N20_terminator_up', 'promoter_N20_terminator_down','seq_altered']]
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['UHA'] + sequencing_primer_df['seq_altered'] + sequencing_primer_df['DHA'] + sequencing_primer_df['promoter_N20_terminator_up'] + sequencing_primer_df['promoter_N20_terminator'] 
        sequencing_primer_df['Region'] = sequencing_primer_df['Name']+';'+sequencing_primer_df['Region']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region')

    return plasmid_sequencing_primer_df

def two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df,promoter_seq,enzyme_df,enzyme_name):

    sgRNA_enzyme_df = enzyme_df[enzyme_df['name']==enzyme_name]
    cut_seq_len = sgRNA_enzyme_df.loc[0,'cut_seq_len']
  

    #N20片段引物退火获得
    temp_sgRNA_df = no_ccdb_uha_dha_sgRNA_df[['Name', 'Region','Target sequence', 'promoter_N20_terminator']]
    temp_sgRNA_df['Region'] = temp_sgRNA_df['Name'] +';'+ temp_sgRNA_df['Region']
    temp_sgRNA_df.drop(columns='Name',inplace=True)
    enzymeCutSeq_and_N20_df = p_d_seq.create_enzymeCutSeq_and_N20(temp_sgRNA_df, promoter_seq, enzyme_cut_len=cut_seq_len)
    #取出重要列
    enzymeCutSeq_and_N20_df = enzymeCutSeq_and_N20_df.rename(columns = {"Region":"ID"})
    enzymeCutSeq_and_N20_df = enzymeCutSeq_and_N20_df[[ "ID","Target sequence","enzymeCutSeq_and_N20"  ]]

    return enzymeCutSeq_and_N20_df

def region_2_distance(sgRNA_plasmid_seq_len, sgRNA_region_json, first_primer_start_position):

    distance_dict = {}
    for k,v in sgRNA_region_json.items():
        arr = v.split(',')
        first = int(arr[0])
        last = int(arr[1])
        if first  > first_primer_start_position:
            distance = first - first_primer_start_position
        else:
            distance = sgRNA_plasmid_seq_len - first_primer_start_position + first
        distance_dict.update({distance:v})
    sorted_distance = sorted(distance_dict.keys())

    last_distance_dict = {}
    for i in range(len(sorted_distance)):
        v = distance_dict[sorted_distance[i]]
        arr = v.split(',')
        first = int(arr[0])
        last = int(arr[1])
        new_len = last - first

        if i != 0:
            v = distance_dict[sorted_distance[i-1]]
            arr = v.split(',')
            first = int(arr[0])
            last = int(arr[1])
            old_len = last - first
            min_distance = sorted_distance[i] - sorted_distance[i-1] -old_len
        else:
            min_distance = sorted_distance[i]  

        max_distance = min_distance + new_len

        last_distance_dict.update({f'primer{i+1}':(min_distance, max_distance)})

    print(last_distance_dict)

    return last_distance_dict

def two_plasmid_system_design_by_user_region(
                                            no_sgRNA_plasmid,
                                            no_ccdb_plasmid,
                                            uha_dha_sgRNA_df, 
                                            enzyme_df,
                                            enzyme_name,
                                            distance_dict,
                                             sgRNA_plasmid_seq,
                                             sgRNA_plasmid_region_seq,
                                             ccdb_plasmid_seq,
                                             ccdb_plasmid_region_seq,
                                             sgRNA_region_json,
                                             ccdb_region_json
                                             ):


    # sgRNA_region_json
    # ccdb_region_json
        #区域转化距离
    # "sgRNA_region_json":{
    #     "region1":"371,570",
    #     "region2":"3572,3770"
    # },
    # distance_dict = {
    #         'distance1':(3000,3200),
    #         'distance2':(3000,3200),
    #     }

    #sgRNA载体质粒                   #启动子终止子若横跨零点有问题
    promoter_terminator_seq = sgRNA_plasmid_region_seq['promoter_seq']+sgRNA_plasmid_region_seq['n20_coordinate_seq']+sgRNA_plasmid_region_seq['terminator_seq']
    promoter_terminator_start = sgRNA_plasmid_seq.find(promoter_terminator_seq)
    promoter_terminator_end = promoter_terminator_start + len(promoter_terminator_seq)

    first_primer_position_in_promoter_terminator = promoter_terminator_seq.find(sgRNA_plasmid_region_seq['terminator_seq'])
    first_primer_start_position = promoter_terminator_start + first_primer_position_in_promoter_terminator
    sgRNA_plasmid_seq_len = len(sgRNA_plasmid_seq)

    sgRNA_distance_dict = region_2_distance(sgRNA_plasmid_seq_len, sgRNA_region_json,first_primer_start_position)
    sgRNA_plasmid_primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(first_primer_start_position, sgRNA_plasmid_seq, sgRNA_distance_dict, 20)
    sgRNA_plasmid_primer_df =  su.result_primer_list_to_df(sgRNA_plasmid_primer_result_list)

    #ccdb载体质粒
    ccdb_start = ccdb_plasmid_seq.find(ccdb_plasmid_region_seq['ccdb'])
    ccdb_end = ccdb_start + len(ccdb_plasmid_region_seq['ccdb'])
    first_primer_start_position = ccdb_end
    ccdb_plasmid_seq_len = len(ccdb_plasmid_seq)
    ccdb_distance_dict = region_2_distance(ccdb_plasmid_seq_len, ccdb_region_json,first_primer_start_position)
    ccdb_plasmid_primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(first_primer_start_position,ccdb_plasmid_seq, ccdb_distance_dict,len(ccdb_plasmid_region_seq['ccdb']))
    ccdb_plasmid_primer_df = su.result_primer_list_to_df(ccdb_plasmid_primer_result_list)

    #sgRNA质粒引物加接头
    sgRNA_plasmid_primer_joint_df = p_d_seq.sgRNA_sgRNAprimer_merge(uha_dha_sgRNA_df, sgRNA_plasmid_primer_df)
    sgRNA_primers_sum=len(sgRNA_plasmid_primer_df)
    sgRNA_plasmid_primer_joint_df = p_d_seq.add_joint_plasmid_primer(enzyme_df,enzyme_name,sgRNA_plasmid_primer_joint_df,sgRNA_primers_sum,primer_type='sgRNA')
    #修改df格式
    sgRNA_plasmid_primer_df = p_d_seq.plasmid_primer(sgRNA_plasmid_primer_joint_df)
    #添加索引
    sgRNA_plasmid_primer_df.insert(0,'index',sgRNA_plasmid_primer_df.index)
    sgRNA_plasmid_primer_df.reset_index(drop=True, inplace=True)
    #给引物加产物
    sgRNA_plasmid_primer_df = p_d_seq.add_product_and_size(no_ccdb_plasmid, sgRNA_plasmid_primer_df, enzyme_df, enzyme_name)

    #ccdB质粒引物加接头
    ccdb_plasmid_primer_joint_df = p_d_seq.sgRNA_sgRNAprimer_merge(uha_dha_sgRNA_df,ccdb_plasmid_primer_df,sgRNA_columns=['Name','Region'])
    ccdb_primers_sum = len(ccdb_plasmid_primer_df)
    ccdb_plasmid_primer_joint_df = p_d_seq.add_joint_plasmid_primer(enzyme_df,enzyme_name,ccdb_plasmid_primer_joint_df,ccdb_primers_sum,primer_type='ccdb')
    ccdb_plasmid_primer_joint_df =  ccdb_plasmid_primer_joint_df.drop(columns='Region').drop_duplicates()
    #修改df格式
    ccdb_plasmid_primer_df = p_d_seq.plasmid_primer(ccdb_plasmid_primer_joint_df)
    #添加索引
    ccdb_plasmid_primer_df.insert(0,'index',ccdb_plasmid_primer_df.index)
    ccdb_plasmid_primer_df.reset_index(drop=True, inplace=True)
    #给引物加产物
    ccdb_plasmid_primer_df = p_d_seq.add_product_and_size(no_sgRNA_plasmid, ccdb_plasmid_primer_df, enzyme_df, enzyme_name='BsaI')

    return  sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df

def two_plasmid_system_design_by_no_user(no_ccdb_uha_dha_sgRNA_df,ccdB_plasmid_backbone,enzyme_df,enzyme_name):
    #设计sgRNA质粒引物
    no_ccdb_primer_template_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','sgRNA_template','Rev Target sequence']]
    no_ccdb_primer_template_df['Region'] = no_ccdb_primer_template_df['Name'] + ';' + no_ccdb_primer_template_df['Region']
    sgRNA_plasmid_primer_df = p_d_seq.design_primer(no_ccdb_primer_template_df,'Region','sgRNA_template','sgRNA')
    no_ccdb_uha_dha_sgRNA_df['Region'] = no_ccdb_uha_dha_sgRNA_df['Name'] +';'+ no_ccdb_uha_dha_sgRNA_df['Region']
    
    #sgRNA质粒引物加接头
    sgRNA_plasmid_primer_df = pd.merge(no_ccdb_uha_dha_sgRNA_df[['Region','Target sequence','Rev Target sequence']],sgRNA_plasmid_primer_df,on=['Region'],how='inner')
    sgRNA_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(sgRNA_plasmid_primer_df,enzyme_df,enzyme_name,stype='sgRNA_plasmid_primer_joint')
    
    #设计ccdB质粒引物
    no_sgRNA_primer_template_df = pd.DataFrame(columns=['plasmid_backbone_primer','plasmid_backbone'],data=[[f'ccdb_plasmid;primer',ccdB_plasmid_backbone]])
    ccdb_plasmid_primer_df = p_d_seq.design_primer(no_sgRNA_primer_template_df,'plasmid_backbone_primer','plasmid_backbone','plasmid')
    
    #ccdb质粒引物加接头
    ccdb_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(ccdb_plasmid_primer_df,enzyme_df,enzyme_name,stype='ccdb_plasmid_primer_joint')
    
    #提取引物的必要部分
    sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
    ccdb_plasmid_p_df = ccdb_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

    return sgRNA_plasmid_p_df, ccdb_plasmid_p_df, sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df

def two_plasmid_system_design_by_user_primer_sgRNA(no_ccdb_plasmid, uha_dha_sgRNA_df, sgRNA_plasmid_backbone, sgRNA_promoter_terminator, sgRNA_primer_json,enzyme_df,enzyme_name,n_20_label,ccdb_label,promoter_terminator_label):
    #定位和检查引物
    primer_position_json, failture_primer = p_d_seq.check_locate_primer(sgRNA_plasmid_backbone, sgRNA_primer_json)

    #对引物进行排序, 确定所有引物
    sgRNA_promoter_terminator_start = sgRNA_plasmid_backbone.find(sgRNA_promoter_terminator)
    sgRNA_primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         sgRNA_primer_json,
                                                         primer_position_json,
                                                         no_ccdb_plasmid,
                                                         sgRNA_plasmid_backbone,
                                                         n_20_label,
                                                         ccdb_label,
                                                         promoter_terminator_label)

    #确定待合并的df
    sgRNA_columns = ['Name','Region','Rev Target sequence']
    temp_sgRNA_df = uha_dha_sgRNA_df[sgRNA_columns]
    temp_sgRNA_df['Region'] = uha_dha_sgRNA_df['Name'] +';'+ uha_dha_sgRNA_df['Region']
    temp_sgRNA_df.drop(columns = 'Name', inplace=True)
    sgRNA_plasmid_primer_joint_df = su.merge_fill_two_df(temp_sgRNA_df, sgRNA_primer_dict_df)
    #给引物加接头
    sgRNA_plasmid_primer_joint_df = p_d_seq.add_joint_plasmid_primer(enzyme_df,enzyme_name,sgRNA_plasmid_primer_joint_df,int(primers_sum/2),primer_type='sgRNA')
    #修改df格式
    sgRNA_plasmid_primer_df = p_d_seq.plasmid_primer(sgRNA_plasmid_primer_joint_df)
    #添加索引
    sgRNA_plasmid_primer_df.insert(0,'index',sgRNA_plasmid_primer_df.index)
    sgRNA_plasmid_primer_df.reset_index(drop=True, inplace=True)   
    #给引物加产物
    sgRNA_plasmid_primer_df = p_d_seq.add_product_and_size(no_ccdb_plasmid, sgRNA_plasmid_primer_df, enzyme_df, enzyme_name)
    return sgRNA_plasmid_primer_df

def two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone, ccdb_primer_json, no_sgRNA_plasmid,enzyme_df,enzyme_name,n_20_label,ccdb_label,promoter_terminator_label):

    #定位和检查引物
    primer_position_json, failture_primer = p_d_seq.check_locate_primer(ccdB_plasmid_backbone, ccdb_primer_json)
    sgRNA_promoter_terminator_start = 0
    #对引物进行排序,确定所有引物
    ccdb_primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         ccdb_primer_json,
                                                         primer_position_json,
                                                         no_sgRNA_plasmid,
                                                         ccdB_plasmid_backbone,
                                                         n_20_label,
                                                         ccdb_label,
                                                         promoter_terminator_label)
    #给引物加接头
    ccdb_plasmid_primer_joint_df = p_d_seq.add_joint_plasmid_primer(enzyme_df,enzyme_name,
                                                                 ccdb_primer_dict_df,
                                                                 int(primers_sum/2),
                                                                 primer_type='ccdb')
    #修改df格式
    ccdb_plasmid_primer_df = p_d_seq.plasmid_primer(ccdb_plasmid_primer_joint_df)
    #添加索引
    ccdb_plasmid_primer_df.insert(0,'index',ccdb_plasmid_primer_df.index)
    ccdb_plasmid_primer_df.reset_index(drop=True, inplace=True)
        
    #给引物加产物
    ccdb_plasmid_primer_df = p_d_seq.add_product_and_size(no_sgRNA_plasmid, ccdb_plasmid_primer_df, enzyme_df,enzyme_name)

    return ccdb_plasmid_primer_df

def two_plasmid_system_design_by_user_primer(uha_dha_sgRNA_df,
                                             sgRNA_primer_json, 
                                             ccdb_primer_json,
                                            sgRNA_plasmid_backbone,
                                            sgRNA_promoter_terminator,
                                            ccdB_plasmid_backbone,
                                            no_sgRNA_plasmid,
                                            enzyme_df,
                                            enzyme_name,
                                            n_20_label,
                                            promoter_terminator_label,
                                            ccdb_label):
    #设计sgRNA primer
    sgRNA_plasmid_primer_df = two_plasmid_system_design_by_user_primer_sgRNA(uha_dha_sgRNA_df, sgRNA_plasmid_backbone, sgRNA_promoter_terminator, sgRNA_primer_json,enzyme_df,enzyme_name,n_20_label,ccdb_label,promoter_terminator_label)

    #设计ccdb primer
    ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone, ccdb_primer_json, no_sgRNA_plasmid,enzyme_df,enzyme_name,n_20_label,ccdb_label,promoter_terminator_label)

    return sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df

def two_plasmid_system_sequencing_design_primer(no_ccdb_uha_dha_sgRNA_df,no_sgRNA_uha_dha_ccdb_df):

    #sgRNA质粒测序
    sgRNA_plasmid_sequencing_primer_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','plasmid','promoter_N20_terminator','promoter_N20_terminator_up','promoter_N20_terminator_down','seq_altered']]
    sgRNA_plasmid_sequencing_primer_df['plasmid_sequencing_region'] = sgRNA_plasmid_sequencing_primer_df['promoter_N20_terminator']
    sgRNA_plasmid_sequencing_primer_df['Region'] = sgRNA_plasmid_sequencing_primer_df['Name']+';'+sgRNA_plasmid_sequencing_primer_df['Region']
    sgRNA_plasmid_sequencing_primer_template = sgRNA_plasmid_sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
    sgRNA_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sgRNA_plasmid_sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region')
    #ccdb质粒载体测序
    ccdb_plasmid_sequencing_primer_df = no_sgRNA_uha_dha_ccdb_df[['Name','Region','UHA','DHA','seq_altered','plasmid']]
    ccdb_plasmid_sequencing_primer_df['plasmid_sequencing_region'] = ccdb_plasmid_sequencing_primer_df['UHA'] + ccdb_plasmid_sequencing_primer_df['seq_altered'] + ccdb_plasmid_sequencing_primer_df['DHA']
    ccdb_plasmid_sequencing_primer_df['Region'] =  ccdb_plasmid_sequencing_primer_df['Name'] + ';' + ccdb_plasmid_sequencing_primer_df['Region']
    ccdb_plasmid_sequencing_primer_template = ccdb_plasmid_sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
    ccdb_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(ccdb_plasmid_sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region')

    return sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df

#genome sequencing primer
def genome_sequencing_design_primer(info_input_df, uha_dha_df):

    #编辑基因组测序
    info_input_df1 = info_input_df[['name','region','seq_uha_max_whole','seq_dha_max_whole','uha_upstream','dha_downstream']].rename(columns={'name':'Name','region':'Region'})
    UHA_DHA_df = pd.merge(info_input_df1,uha_dha_df)
    UHA_DHA_df['sequencing_template'] =UHA_DHA_df['uha_upstream'] + UHA_DHA_df['seq_uha_max_whole']+UHA_DHA_df['seq_altered']+UHA_DHA_df['seq_dha_max_whole'] + UHA_DHA_df['dha_downstream'] 
    UHA_DHA_df['sequencing_region'] = UHA_DHA_df['UHA']+UHA_DHA_df['seq_altered']+UHA_DHA_df['DHA']
    UHA_DHA_df['Region'] =  UHA_DHA_df['Name']+';'+UHA_DHA_df['Region']
    genome_sequencing_primer_df = p_d_seq.create_sequencing_primer(UHA_DHA_df,sr,'sequencing_template','sequencing_region')

    return genome_sequencing_primer_df

# 
def execute_one_plasmid_system(gb_path, info_df, info_input_df, uha_dha_df, uha_dha_sgRNA_df, uha_dha_info_primer_df, uha_dha_primer_df,enzyme_df,enzyme_name,output):
    gb = SeqIO.read(gb_path, "genbank")
    gb_name = gb.name

    #设计pcr引物
    uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_primer_df,type_kind = one_plasmid_system_pcr_design_primer(gb_path,
                                                                                                                                    info_df,
                                                                                                                                    uha_dha_sgRNA_df,
                                                                                                                                    uha_dha_info_primer_df,
                                                                                                                                    uha_dha_primer_df,
                                                                                                                                    enzyme_df,enzyme_name)
    #设计质粒测序引物
    plasmid_sequencing_primer_df = one_plasmid_system_sequencing_design_primer(type_kind,uha_dha_sgRNA_df)
    #设计基因组测序引物
    genome_sequencing_primer_df = genome_sequencing_design_primer(info_input_df, uha_dha_df)
    #标准化重命名,uha_primer_df,dha_primer_df,sgRNA_primer_df,plasmid_sequencing_primer_df,genome_sequencing_primer_df
    uha_primer_df,dha_primer_df,sgRNA_primer_df,plasmid_sequencing_primer_df,genome_sequencing_primer_df,plasmid_backbone_p_df,seq_altered_p_df = su.normal_primer_rename(
                                                                                                                                                        uha_primer_df,
                                                                                                                                                        dha_primer_df,
                                                                                                                                                        n20up_primer_p_df,
                                                                                                                                                        plasmid_sequencing_primer_df,                                                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                        genome_sequencing_primer_df,
                                                                                                                                                        n20down_primer_p_df,                                                                                                                                         
                                                                                                                                                        seq_altered_primer_df)
    #输出引物  
    with pd.ExcelWriter(output+'one_plasmid_design_result.xlsx') as writer:  
            uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            sgRNA_primer_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')  
            plasmid_backbone_p_df.to_excel(writer,sheet_name = 'Primer_plasmid',index_label='No.')
            seq_altered_p_df.to_excel(writer,sheet_name = 'Seq_altered',index_label='No.')
            plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P',index_label='No.')
            genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.') 

    return output + 'one_plasmid_design_result.xlsx'                                                                                                                                       

def execute_two_plasmid_system(
                                 info_df,
                                 uha_dha_info_primer_df,  
                                 uha_dha_df,
                                 info_input_df,
                                 no_ccdb_plasmid,
                                 no_sgRNA_plasmid,
                                 uha_dha_sgRNA_df,
                                 plasmid_primer_desgin_type,
                                 enzyme_df,
                                 enzyme_name,
                                 sgRNA_primer_json,
                                 ccdb_primer_json,
                                 sgRNA_region_json,
                                 ccdb_region_json,
                                 ccdb_label,
                                 promoter_terminator_label,
                                 n_20_label,
                                 output
                                 ):
    no_ccdb_uha_dha_sgRNA_df, sgRNA_plasmid_backbone, promoter_seq, terminator_seq, sgRNA_promoter_terminator = p_d_seq.create_new_plasmid(no_ccdb_plasmid,uha_dha_sgRNA_df.copy(),ccdb_label, promoter_terminator_label, n_20_label)
    no_sgRNA_uha_dha_ccdb_df, ccdB_plasmid_backbone, ccdB_promoter_terminator_up_seq = p_d_seq.create_new_plasmid(no_sgRNA_plasmid,uha_dha_sgRNA_df.copy(),ccdb_label, promoter_terminator_label, n_20_label)
    #酶切退火方式
    enzymeCutSeq_and_N20_df = two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df,promoter_seq,enzyme_df,enzyme_name)

    #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物
    if  plasmid_primer_desgin_type == 1:
        #对质粒进行分割
        #初步分割
        sgRNA_plasmid_seq, sgRNA_plasmid_region_seq = p_d_seq.plasmid_region_division_by_labels(no_ccdb_plasmid, ccdb_label, promoter_terminator_label, n_20_label)
        ccdb_plasmid_seq, ccdb_plasmid_region_seq = p_d_seq.plasmid_region_division_by_labels(no_sgRNA_plasmid, ccdb_label, promoter_terminator_label, n_20_label)

        #第一段距离：sgRNA终止子距离第一段区域的最小、最大距离
        #第二段距离：是第一段区域距离第二段区域的最小、最大距离
        #。。。。。。。。。
        #最后一段距离:最后一个区域距离sgRNA启动子距离

        #区域转换距离

        distance_dict = {
            'distance1':(3000,3200),
            'distance2':(3000,3200),
        }
        #根据区域设计引物
        sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_region(
                                                                                    no_sgRNA_plasmid,
                                                                                    no_ccdb_plasmid,
                                                                                    uha_dha_sgRNA_df,
                                                                                    enzyme_df,
                                                                                    enzyme_name,
                                                                                    distance_dict, 
                                                                                    sgRNA_plasmid_seq,
                                                                                    sgRNA_plasmid_region_seq,
                                                                                    ccdb_plasmid_seq,
                                                                                    ccdb_plasmid_region_seq,
                                                                                    sgRNA_region_json,
                                                                                    ccdb_region_json)

    elif    plasmid_primer_desgin_type == 2:
        sgRNA_plasmid_p_df, ccdb_plasmid_p_df, sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df = two_plasmid_system_design_by_no_user(no_ccdb_uha_dha_sgRNA_df, ccdB_plasmid_backbone,enzyme_df,enzyme_name)
    elif    plasmid_primer_desgin_type == 3:
         #用户提供正义链上的引物
        # 引物 AACTATTTATCCAGTTGGTACAAAC
        # sgRNA_primer_json = {
        #     'primer3':'AACTATTTATCCAGTTGGTACAAAC',
        #     'primer5':'GAAGATAACGAACAAAAACAATTGT'
        # }
        # ccdb_primer_json = {
        #     'primer3':"AACTGATTCAGTCTGATTTCGCGGT",
        #     'primer5':"CCCTCTAATCGAAACTAATGGGGA"
        # }
        #判断引物类型，设计引物
        #sgRNA_primer
        sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer(uha_dha_sgRNA_df,
                                                                                                    sgRNA_primer_json, 
                                                                                                    ccdb_primer_json,
                                                                                                    sgRNA_plasmid_backbone,
                                                                                                    sgRNA_promoter_terminator,
                                                                                                    ccdB_plasmid_backbone,
                                                                                                    no_sgRNA_plasmid,
                                                                                                    enzyme_df,
                                                                                                    enzyme_name,
                                                                                                    n_20_label,
                                                                                                    promoter_terminator_label,
                                                                                                    ccdb_label)
    
    #设计seq_altered_primer，seq_altered > 120
    seq_altered_primer_template =  info_df[info_df.seq_altered.apply(lambda x:len(x)>120)][['Name','Region','seq_altered']]
    seq_altered_primer_template['Region'] = seq_altered_primer_template['Name'] + seq_altered_primer_template['Region']
    seq_altered_primer_df = p_d_seq.design_primer(seq_altered_primer_template,'Region','seq_altered',stype='seq_altered')

    #给引物加接头
    #给同源臂引物加接头
    uha_dha_primer_df = p_d_seq.add_joint_sgRNA_primer(uha_dha_info_primer_df, enzyme_df, enzyme_name, ccdB_plasmid_backbone, stype = 'u_d_primer_joint')
    #seq_altered_primer加接头
    seq_altered_primer_df = p_d_seq.add_joint_sgRNA_primer(seq_altered_primer_df, enzyme_df, enzyme_name, stype='seq_altered_primer_joint')

    #分别提取，修饰后的uha、dha引物
    in_col = ['Name','Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']
    ou_col = ['Name','Region',"u_primer_f_seq_(5'-3')","u_primer_r_seq_(5'-3')",'UHA','UHA_size',"d_primer_f_seq_(5'-3')","d_primer_r_seq_(5'-3')",'DHA','DHA_size']
    uha_primer_df, dha_primer_df = p_d_seq.create_uha_dha_primer_df(uha_dha_primer_df, in_col, ou_col)
    #提取变化序列引物
    seq_altered_p_df = seq_altered_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
    if 'index' in sgRNA_plasmid_primer_df.columns and 'index' in ccdb_plasmid_primer_df.columns:
        sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df
        ccdb_plasmid_p_df = ccdb_plasmid_primer_df
    else:
        sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
        ccdb_plasmid_p_df = ccdb_plasmid_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

    #设计质粒测序引物
    sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df = two_plasmid_system_sequencing_design_primer(no_ccdb_uha_dha_sgRNA_df,no_sgRNA_uha_dha_ccdb_df)
    #设计基因组测序引物
    genome_sequencing_primer_df = genome_sequencing_design_primer(info_input_df, uha_dha_df)

    #标准化，重命名
    df_common_list = su.rename_common_primer_df(sgRNA_plasmid_p_df,ccdb_plasmid_p_df,seq_altered_p_df)
    sgRNA_plasmid_p_df, ccdb_plasmid_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1], df_common_list[2]
    uha_primer_df, dha_primer_df = su.rename_u_d_primer_df(uha_primer_df,dha_primer_df)
    df_sequencing_list = su.rename_sequencing_primer_df(sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df,genome_sequencing_primer_df)
    sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df, genome_sequencing_primer_df = df_sequencing_list[0], df_sequencing_list[1], df_sequencing_list[2]

    #输出引物  
    with pd.ExcelWriter(output+'two_plasmid_design_result.xlsx') as writer:  
        uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
        dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
        sgRNA_plasmid_p_df.to_excel(writer,sheet_name = 'sgRNA_plasmid_fragment',index_label='No.')  
        ccdb_plasmid_p_df.to_excel(writer,sheet_name = 'uha_dha_plasmid',index_label='No.')
        seq_altered_p_df.to_excel(writer,sheet_name = 'Seq_altered',index_label='No.')
        sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
        ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
        genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')

    return output + 'two_plasmid_design_result.xlsx'


def main(data):

    chopchop_input = data['chopchop_input']
    # enzyme_path = parent_base_path +'/'+ data['enzyme_path']
    sgRNA_result_path = data['sgRNA_result_path']
    one_plasmid_file_path = data['one_plasmid_file_path']
    no_ccdb_plasmid = data['no_ccdb_plasmid']
    no_sgRNA_plasmid = data['no_sgRNA_plasmid']
    
    #plasmid label
    ccdb_label = data['plasmid_label']['ccdb_label']
    promoter_terminator_label = data['plasmid_label']['promoter_terminator_label']
    n_20_label = data['plasmid_label']['n_20_label']   
 
    #primer
    sgRNA_primer_json = data['sgRNA_primer_json']
    ccdb_primer_json = data['ccdb_primer_json']

    sgRNA_region_json = data['sgRNA_region_json']
    ccdb_region_json = data['ccdb_region_json']

    #配置引物参数
    config.S_GLOBAL_ARGS = data['S_GLOBAL_ARGS']
    config.Q_ARGS = data['Q_ARGS']


    #配置输出参数
    output = data['edit_sequence_design_workdir']
    if not os.path.exists(output):
        os.makedirs(output)
    # 1.read 编辑序列信息
    info_input_df = pd.read_csv(chopchop_input)    
    # 2.#读取用户填的酶
    # enzyme={
    #     "enzyme_name":"BsaI",  
    #     "protective_base":"CCA",
    #     "recognition_seq":"GGTCTC",
    #     "cut_seq_len":4,
    #     "gap_len":1    
    # }
    enzyme = data['enzyme']
    enzyme_df = pd.DataFrame([enzyme])
    enzyme_name =data['enzyme']['name']
    # enzyme_df = pd.read_csv(enzyme_path) 

    # 3.提取用户选择的sgRNA
    sgRNA = p_d_seq.extract_sgRNA_from_chopchop(sgRNA_result_path)

    # 4.设计源生同源臂引物
    uha_dha_primer_df = extract_uha_dha_primer(info_input_df, sgRNA)

    # 5.提取同源臂
    uha_dha_info_primer_df, uha_dha_df, uha_dha_sgRNA_df, info_df = extract_uha_dha(info_input_df,uha_dha_primer_df,sgRNA)




    if  one_plasmid_file_path != '' and no_ccdb_plasmid == '' and no_sgRNA_plasmid == '':
        plasmid_system_type =1
    elif    one_plasmid_file_path =='' and no_ccdb_plasmid != '' and no_sgRNA_plasmid != '':
        plasmid_system_type = 2
    else:
        plasmid_system_type = 0
    
    print('---------------------------------------------------------',plasmid_system_type) 

    if plasmid_system_type == 1:
        # 6.执行单质粒系统
        one_plasmid_output_path = execute_one_plasmid_system(one_plasmid_file_path, info_df, info_input_df, uha_dha_df, uha_dha_sgRNA_df, uha_dha_info_primer_df, uha_dha_primer_df,enzyme_df,enzyme_name,output)
        return one_plasmid_output_path  
    elif plasmid_system_type == 2:
        # 7.执行双质粒系统
        #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物
       
        plasmid_primer_desgin_type = 2
        if sgRNA_primer_json !={} and ccdb_primer_json !={} and sgRNA_region_json =={} and ccdb_region_json =={}:
            plasmid_primer_desgin_type = 3
        elif sgRNA_region_json !={} and ccdb_region_json !={} and sgRNA_primer_json =={} and ccdb_primer_json =={}:
            plasmid_primer_desgin_type = 1
        else:
            plasmid_primer_desgin_type = 2 
        print('-----------------------------------------',plasmid_primer_desgin_type)    

        two_plasmid_output_path = execute_two_plasmid_system(
                                    info_df,
                                    uha_dha_info_primer_df,
                                    uha_dha_df,
                                    info_input_df,
                                    no_ccdb_plasmid,
                                    no_sgRNA_plasmid,
                                    uha_dha_sgRNA_df,
                                    plasmid_primer_desgin_type,
                                    enzyme_df,
                                    enzyme_name,
                                    sgRNA_primer_json,
                                    ccdb_primer_json,
                                    sgRNA_region_json,
                                    ccdb_region_json,
                                    ccdb_label,
                                    promoter_terminator_label,
                                    n_20_label,
                                    output
                                    )
        return two_plasmid_output_path
    elif plasmid_system_type ==0:

        one_plasmid_output_path = execute_one_plasmid_system(one_plasmid_file_path, info_df, info_input_df, uha_dha_df, uha_dha_sgRNA_df, uha_dha_info_primer_df, uha_dha_primer_df,enzyme_df,enzyme_name,output)

        plasmid_primer_desgin_type = 2
        if sgRNA_primer_json !={} and ccdb_primer_json !={} and sgRNA_region_json =={} and ccdb_region_json =={}:
            plasmid_primer_desgin_type = 3
        elif sgRNA_region_json !={} and ccdb_region_json !={} and sgRNA_primer_json =={} and ccdb_primer_json =={}:
            plasmid_primer_desgin_type = 1
        else:
            plasmid_primer_desgin_type = 2 
        print('-----------------------------------------',plasmid_primer_desgin_type)    

        two_plasmid_output_path = execute_two_plasmid_system(
                                    info_df,
                                    uha_dha_info_primer_df,
                                    uha_dha_df,
                                    info_input_df,
                                    no_ccdb_plasmid,
                                    no_sgRNA_plasmid,
                                    uha_dha_sgRNA_df,
                                    plasmid_primer_desgin_type,
                                    enzyme_df,
                                    enzyme_name,
                                    sgRNA_primer_json,
                                    ccdb_primer_json,
                                    sgRNA_region_json,
                                    ccdb_region_json,
                                    ccdb_label,
                                    promoter_terminator_label,
                                    n_20_label,
                                    output
                                    )
    return one_plasmid_output_path,two_plasmid_output_path
     
if __name__ == '__main__':

    #read json
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help='input params file', required=True) 
    args = parser.parse_args()
    input_path =  args.input
    with open(input_path, "r") as f:
        data = json.load(f)

    main(data)

    
