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

from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
# from loguru import logger

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

def one_plasmid_system_design_by_user_primer(primer_json,primer_position_json,gb_path,plasmid_backbone,enzyme_df,enzyme_name):

    sgRNA_promoter_terminator_start = 0
    #对引物进行排序,确定所有引物
    primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         primer_json,
                                                         primer_position_json,
                                                         gb_path,
                                                         plasmid_backbone,
                                                         plasmid_type='one')   
    #给引物加接头
    plasmid_primer_joint_df = p_d_seq.add_joint_plasmid_primer(enzyme_df,enzyme_name,
                                                                 primer_dict_df,
                                                                 int(primers_sum/2),
                                                                 primer_type='ccdb')
    #修改df格式
    plasmid_primer_df = p_d_seq.plasmid_primer(plasmid_primer_joint_df)
    #添加索引
    plasmid_primer_df.insert(0,'index',plasmid_primer_df.index)
    plasmid_primer_df.reset_index(drop=True, inplace=True)
        
    #给引物加产物
    plasmid_primer_df = p_d_seq.add_product_and_size(gb_path, plasmid_primer_df, enzyme_df,enzyme_name,seq=plasmid_backbone)

    return  plasmid_primer_df

def one_plasmid_system_design_by_user_region(uha_dha_sgRNA_df,sgRNA_region_seq_json,one_plasmid_file_path,sgRNA_plasmid_backbone,enzyme_df,enzyme_name):
    
    #序列转换成坐标
    region_json = su.convert_seq_cor(one_plasmid_file_path, sgRNA_region_seq_json,seq=sgRNA_plasmid_backbone)
    #坐标转转换成距离
    distance_dict = region_2_distance(len(sgRNA_plasmid_backbone), region_json, 0)

    #根据区域设计引物
    primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(0,sgRNA_plasmid_backbone,distance_dict)
    plasmid_primer_df = su.result_primer_list_to_df(primer_result_list)

    #加接头
    plasmid_primer_joint_df = p_d_seq.sgRNA_sgRNAprimer_merge(uha_dha_sgRNA_df,plasmid_primer_df,sgRNA_columns=['Name','Region'])
    primers_sum = len(plasmid_primer_df)
    plasmid_primer_joint_df = p_d_seq.add_joint_plasmid_primer(enzyme_df, enzyme_name, plasmid_primer_joint_df, primers_sum, primer_type='ccdb')
    plasmid_primer_joint_df = plasmid_primer_joint_df.drop(columns='Region').drop_duplicates()
    #修改df格式
    plasmid_primer_df = p_d_seq.plasmid_primer(plasmid_primer_joint_df)
    #添加索引
    plasmid_primer_df.insert(0,'index',plasmid_primer_df.index)
    plasmid_primer_df.reset_index(drop=True, inplace=True)

    plasmid_primer_df = p_d_seq.add_product_and_size(one_plasmid_file_path, plasmid_primer_df, enzyme_df, enzyme_name=enzyme_name,seq=sgRNA_plasmid_backbone)

    return  plasmid_primer_df

#one plasmid pcr primer     
def one_plasmid_system_pcr_design_primer(gb_path,
                                         info_df,
                                         uha_dha_sgRNA_df,
                                         uha_dha_info_primer_df,
                                         uha_dha_primer_df,
                                         enzyme_df,
                                         enzyme_name,
                                         plasmid_primer_desgin_type,
                                         seq_json):
    
    #创建新的质粒   
    uha_dha_sgRNA_df, promoter_terminator_up_promoter_seq, promoter_terminator_down_terminator_seq, type_kind  = p_d_seq.create_new_plasmid(gb_path, uha_dha_sgRNA_df.copy(), ccdb_label='ccdB', promoter_terminator_label='gRNA', n_20_label='N20')

    #设计sgRNA、ccdb质粒引物
    n20up_primer_template = uha_dha_sgRNA_df[['Name','Region','n20_up_template','Target sequence','Rev Target sequence']]
    n20up_primer_template['Region'] = n20up_primer_template['Name'] +';'+ n20up_primer_template['Region']
    n20up_primer_df = p_d_seq.design_primer(n20up_primer_template,'Region','n20_up_template','sgRNA')
    n20up_primer_df = pd.merge(n20up_primer_template[['Region','Target sequence','Rev Target sequence']],n20up_primer_df,on=['Region'],how='inner')

    #ccdb、sgrna质粒引物加接头
    n20up_primer_df = p_d_seq.add_joint_sgRNA_primer(n20up_primer_df,enzyme_df, enzyme_name,'',stype='n20up_primer_joint')


     #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物  
    if plasmid_primer_desgin_type == 2: 
        n20down_primer_template = uha_dha_sgRNA_df[['Name','Region','n20_down_template','Target sequence','Rev Target sequence']]
        n20down_primer_template['Region'] = n20down_primer_template['Name'] +';'+ n20down_primer_template['Region']
        n20down_primer_df = p_d_seq.design_primer(n20down_primer_template,'Region','n20_down_template','sgRNA')
        n20down_primer_df = pd.merge(n20down_primer_template[['Region','Target sequence','Rev Target sequence']],n20down_primer_df,on=['Region'],how='inner')
        #加接头
        n20down_primer_df = p_d_seq.add_joint_sgRNA_primer(n20down_primer_df,enzyme_df,enzyme_name,'',stype='n20down_primer_joint')

    elif plasmid_primer_desgin_type == 1:
        plasmid_backbone = promoter_terminator_down_terminator_seq
        n20down_primer_df =  one_plasmid_system_design_by_user_region(uha_dha_sgRNA_df,seq_json,gb_path,plasmid_backbone,enzyme_df,enzyme_name)
    elif plasmid_primer_desgin_type == 3:
        plasmid_backbone = promoter_terminator_down_terminator_seq
        primer_position_json, sgRNA_failture_primer = p_d_seq.check_locate_primer(plasmid_backbone, seq_json)
        n20down_primer_df = one_plasmid_system_design_by_user_primer(seq_json,primer_position_json,gb_path,plasmid_backbone,enzyme_df,enzyme_name)


    seq_altered_primer_template = info_df[info_df.seq_altered.apply(lambda x:len(x)>120)][['Name','Region','seq_altered']]
    seq_altered_primer_template['Region'] = seq_altered_primer_template['Name'] +';'+ seq_altered_primer_template['Region']
    seq_altered_primer_df = p_d_seq.design_primer(seq_altered_primer_template,'Region','seq_altered',stype='seq_altered')

    #给同源臂引物加接头:uha取promoter_terminator_down_terminator_seq尾部反义4bp，dha取头promoter_terminator_up_promoter_seq正义4bp
    uha_dha_primer_df = p_d_seq.add_joint_sgRNA_primer(uha_dha_info_primer_df,enzyme_df,enzyme_name,promoter_terminator_down_terminator_seq, promoter_terminator_up_promoter_seq, stype = 'u_d_primer_joint')

    #seq_altered_primer加接头
    seq_altered_primer_df = p_d_seq.add_joint_sgRNA_primer(seq_altered_primer_df,enzyme_df,enzyme_name,'',stype='seq_altered_primer_joint')

    #分别提取，修饰后的uha、dha引物
    in_col = ['Name','Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']
    ou_col = ['Name','Region',"u_primer_f_seq_(5'-3')","u_primer_r_seq_(5'-3')",'UHA','UHA_size',"d_primer_f_seq_(5'-3')","d_primer_r_seq_(5'-3')",'DHA','DHA_size']
    uha_primer_df, dha_primer_df = p_d_seq.create_uha_dha_primer_df(uha_dha_primer_df, in_col, ou_col)

    n20down_primer_p_df = n20down_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]
    n20up_primer_p_df = n20up_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]
    seq_altered_primer_df = seq_altered_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]

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
        plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')

    elif type_kind == 1:
        #载体测序
        #uha_dha测序
        sequencing_primer_df=uha_dha_sgRNA_df[['Name', 'Region', 'UHA', 'UHA_size', 'DHA', 'DHA_size',
                                    'plasmid', 'promoter_N20_terminator',
                                    'promoter_N20_terminator_up', 'promoter_N20_terminator_down','seq_altered']]
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['UHA'] + sequencing_primer_df['seq_altered'] + sequencing_primer_df['DHA']
        sequencing_primer_df['Region'] = sequencing_primer_df['Name']+';'+sequencing_primer_df['Region']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df1 = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')
        #promoter_N20_terminator测序
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['promoter_N20_terminator']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df2 = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')
        plasmid_sequencing_primer_df = p_d_seq.merge_sequencing_result(plasmid_sequencing_primer_df1, plasmid_sequencing_primer_df2)
        plasmid_sequencing_primer_df2.to_csv('plasmid_sequencing_primer_df2.csv',index=False)
    elif type_kind == 2:  
        #载体测序  
        #uha_dha测序
        sequencing_primer_df=uha_dha_sgRNA_df[[ 'Name', 'Region', 'UHA', 'UHA_size', 'DHA', 'DHA_size',
                                                'plasmid', 'promoter_N20_terminator',
                                                'promoter_N20_terminator_up', 'promoter_N20_terminator_down','seq_altered']]
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['UHA'] + sequencing_primer_df['seq_altered'] + sequencing_primer_df['DHA'] + sequencing_primer_df['promoter_N20_terminator_up'] + sequencing_primer_df['promoter_N20_terminator'] 
        sequencing_primer_df['Region'] = sequencing_primer_df['Name']+';'+sequencing_primer_df['Region']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')

    temp = uha_dha_sgRNA_df[['Name','Region','UHA','DHA','promoter_N20_terminator','seq_altered']]
    temp['Region'] = temp['Name']+';'+temp['Region']
    temp = temp.drop(columns='Name')
    sequencing_primer_template = pd.merge(sequencing_primer_template,temp,on=['Region'],how='inner')

    return plasmid_sequencing_primer_df, sequencing_primer_template  

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

def compute_distance_in_gb(gb_len,first,last):
    if last < first:
        seq_len = gb_len - first + last
    else:
        seq_len = last - first
    return seq_len

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
        new_len = compute_distance_in_gb(sgRNA_plasmid_seq_len,first,last)
        
        if i != 0:
            v = distance_dict[sorted_distance[i-1]]
            arr = v.split(',')
            first = int(arr[0])
            last = int(arr[1])
            old_len = compute_distance_in_gb(sgRNA_plasmid_seq_len,first,last)
                
            min_distance = sorted_distance[i] - sorted_distance[i-1] -old_len
        else:
            min_distance = sorted_distance[i]  

        max_distance = min_distance + new_len

        last_distance_dict.update({f'primer{i+1}':(min_distance, max_distance)})

    print(last_distance_dict)

    return last_distance_dict

def two_plasmid_system_design_by_user_region(
                                            no_ccdb_uha_dha_sgRNA_df,
                                            ccdB_plasmid_backbone,
                                            no_sgRNA_plasmid,
                                            no_ccdb_plasmid,
                                            uha_dha_sgRNA_df, 
                                            enzyme_df,
                                            enzyme_name,
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

    if sgRNA_region_json != {}:
        sgRNA_plasmid_primer_df = two_plasmid_system_design_by_user_sgRNA_region(uha_dha_sgRNA_df,sgRNA_plasmid_seq,sgRNA_plasmid_region_seq,sgRNA_region_json,no_ccdb_plasmid,enzyme_df,enzyme_name)
    else:
        sgRNA_plasmid_primer_df = two_plasmid_system_design_by_sgRNA_no_user(no_ccdb_uha_dha_sgRNA_df, uha_dha_sgRNA_df, enzyme_df, enzyme_name)

    if ccdb_region_json != {}:
        ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_ccdb_region(uha_dha_sgRNA_df,no_sgRNA_plasmid,ccdb_plasmid_seq,ccdb_plasmid_region_seq,ccdb_region_json, enzyme_df, enzyme_name)
    else:
        ccdb_plasmid_primer_df = two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name)

    return  sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df   


def two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name):
        
        #执行用户什么都不指定，设计ccdb质粒引物   
        ccdb_primer_template_df = pd.DataFrame(columns=['plasmid_backbone_primer','plasmid_backbone'],data=[[f'ccdb_plasmid;primer',ccdB_plasmid_backbone]])
        ccdb_plasmid_primer_df = p_d_seq.design_primer(ccdb_primer_template_df,'plasmid_backbone_primer','plasmid_backbone','plasmid')
        #ccdb质粒引物加接头
        ccdb_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(ccdb_plasmid_primer_df,enzyme_df,enzyme_name,stype='ccdb_plasmid_primer_joint')
        #提取引物的必要部分
        ccdb_plasmid_primer_df = ccdb_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]    

        return ccdb_plasmid_primer_df

def two_plasmid_system_design_by_sgRNA_no_user(no_ccdb_uha_dha_sgRNA_df, uha_dha_sgRNA_df, enzyme_df, enzyme_name):
        
        #执行用户什么都不指定，设计sgRNA质粒引物
        sgRNA_primer_template_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','sgRNA_template','Rev Target sequence']]
        sgRNA_primer_template_df['Region'] = sgRNA_primer_template_df['Name'] +';'+ sgRNA_primer_template_df['Region']
        sgRNA_plasmid_primer_df = p_d_seq.design_primer(sgRNA_primer_template_df,'Region','sgRNA_template','sgRNA')
        uha_dha_sgRNA_df['Region'] = uha_dha_sgRNA_df['Name'] + ';' + uha_dha_sgRNA_df['Region']
        #sgRNA质粒引物加接头
        sgRNA_plasmid_primer_df = pd.merge(uha_dha_sgRNA_df[['Region','Target sequence','Rev Target sequence']],sgRNA_plasmid_primer_df,on=['Region'],how='inner')
        sgRNA_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(sgRNA_plasmid_primer_df,enzyme_df,enzyme_name,stype='sgRNA_plasmid_primer_joint')
        #提取关键信息
        sgRNA_plasmid_primer_df = sgRNA_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

        return sgRNA_plasmid_primer_df


def two_plasmid_system_design_by_user_sgRNA_region(uha_dha_sgRNA_df,sgRNA_plasmid_seq,sgRNA_plasmid_region_seq,sgRNA_region_json,no_ccdb_plasmid,enzyme_df,enzyme_name):
        
        #sgRNA载体质粒                   #启动子终止子若横跨零点有问题
        promoter_terminator_seq = sgRNA_plasmid_region_seq['promoter_seq'] + sgRNA_plasmid_region_seq['n20_coordinate_seq'] + sgRNA_plasmid_region_seq['terminator_seq']
        promoter_terminator_start = sgRNA_plasmid_seq.find(promoter_terminator_seq)
        promoter_terminator_end = promoter_terminator_start + len(promoter_terminator_seq)

        first_primer_position_in_promoter_terminator = promoter_terminator_seq.find(sgRNA_plasmid_region_seq['terminator_seq'])
        first_primer_start_position = promoter_terminator_start + first_primer_position_in_promoter_terminator
        sgRNA_plasmid_seq_len = len(sgRNA_plasmid_seq)

        sgRNA_distance_dict = region_2_distance(sgRNA_plasmid_seq_len, sgRNA_region_json,first_primer_start_position)
        sgRNA_plasmid_primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(first_primer_start_position, sgRNA_plasmid_seq, sgRNA_distance_dict, 20)
        sgRNA_plasmid_primer_df =  su.result_primer_list_to_df(sgRNA_plasmid_primer_result_list)

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

        return sgRNA_plasmid_primer_df


def two_plasmid_system_design_by_user_ccdb_region(uha_dha_sgRNA_df,no_sgRNA_plasmid,ccdb_plasmid_seq,ccdb_plasmid_region_seq,ccdb_region_json, enzyme_df, enzyme_name):

        #ccdb载体质粒
        ccdb_start = ccdb_plasmid_seq.find(ccdb_plasmid_region_seq['ccdb'])
        ccdb_end = ccdb_start + len(ccdb_plasmid_region_seq['ccdb'])
        first_primer_start_position = ccdb_end
        ccdb_plasmid_seq_len = len(ccdb_plasmid_seq)
        ccdb_distance_dict = region_2_distance(ccdb_plasmid_seq_len, ccdb_region_json,first_primer_start_position)
        ccdb_plasmid_primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(first_primer_start_position,ccdb_plasmid_seq, ccdb_distance_dict,len(ccdb_plasmid_region_seq['ccdb']))
        ccdb_plasmid_primer_df = su.result_primer_list_to_df(ccdb_plasmid_primer_result_list)
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
        ccdb_plasmid_primer_df = p_d_seq.add_product_and_size(no_sgRNA_plasmid, ccdb_plasmid_primer_df, enzyme_df, enzyme_name=enzyme_name)

        return ccdb_plasmid_primer_df


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

def two_plasmid_system_design_by_user_primer_sgRNA( no_ccdb_plasmid,
                                                    uha_dha_sgRNA_df,
                                                    sgRNA_plasmid_backbone,
                                                    sgRNA_promoter_terminator,
                                                    sgRNA_primer_json,
                                                    sgRNA_primer_position_json,
                                                    enzyme_df,
                                                    enzyme_name,
                                                    n_20_label,
                                                    ccdb_label,
                                                    promoter_terminator_label):


    #对引物进行排序, 确定所有引物
    sgRNA_promoter_terminator_start = sgRNA_plasmid_backbone.find(sgRNA_promoter_terminator)
    sgRNA_primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         sgRNA_primer_json,
                                                         sgRNA_primer_position_json,
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

def two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone, ccdb_primer_json, ccdb_primer_position_json, no_sgRNA_plasmid,enzyme_df,enzyme_name,n_20_label,ccdb_label,promoter_terminator_label):

    
    sgRNA_promoter_terminator_start = 0
    #对引物进行排序,确定所有引物
    ccdb_primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         ccdb_primer_json,
                                                         ccdb_primer_position_json,
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

def two_plasmid_system_design_by_user_primer(
                                            no_ccdb_uha_dha_sgRNA_df,
                                            uha_dha_sgRNA_df,
                                             sgRNA_primer_json, 
                                             ccdb_primer_json,
                                            sgRNA_plasmid_backbone,
                                            sgRNA_promoter_terminator,
                                            ccdB_plasmid_backbone,
                                            no_ccdb_plasmid,
                                            no_sgRNA_plasmid,
                                            enzyme_df,
                                            enzyme_name,
                                            n_20_label,
                                            promoter_terminator_label,
                                            ccdb_label):

    #定位和检查引物
    sgRNA_primer_position_json = {}
    ccdb_primer_position_json = {}

    if sgRNA_primer_json !={}:
        sgRNA_primer_position_json, sgRNA_failture_primer = p_d_seq.check_locate_primer(sgRNA_plasmid_backbone, sgRNA_primer_json)
    else:
        sgRNA_plasmid_primer_df = two_plasmid_system_design_by_sgRNA_no_user(no_ccdb_uha_dha_sgRNA_df, uha_dha_sgRNA_df, enzyme_df, enzyme_name)
        
    if ccdb_primer_json !={}:  
        ccdb_primer_position_json, ccdb_failture_primer = p_d_seq.check_locate_primer(ccdB_plasmid_backbone, ccdb_primer_json)
    else:
        ccdb_plasmid_primer_df = two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name)



    if  sgRNA_primer_position_json == {} and sgRNA_primer_json !={}:
        return sgRNA_failture_primer
    elif sgRNA_primer_position_json !={} and sgRNA_primer_json != {}:
        #设计sgRNA primer   
        sgRNA_plasmid_primer_df = two_plasmid_system_design_by_user_primer_sgRNA( 
                                                                                no_ccdb_plasmid,
                                                                                uha_dha_sgRNA_df,
                                                                                sgRNA_plasmid_backbone,
                                                                                sgRNA_promoter_terminator,
                                                                                sgRNA_primer_json,
                                                                                sgRNA_primer_position_json,
                                                                                enzyme_df,
                                                                                enzyme_name,
                                                                                n_20_label,
                                                                                ccdb_label,
                                                                                promoter_terminator_label)
    
    if  ccdb_primer_position_json == {} and ccdb_primer_json !={}:
        return ccdb_failture_primer
    elif ccdb_primer_position_json != {} and ccdb_primer_json != {}:
        #设计ccdb primer
        ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone,
                                                                                ccdb_primer_json,
                                                                                ccdb_primer_position_json,
                                                                                no_sgRNA_plasmid,
                                                                                enzyme_df,
                                                                                enzyme_name,
                                                                                n_20_label,
                                                                                ccdb_label,
                                                                                promoter_terminator_label)
              
    return sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df

def two_plasmid_system_sequencing_design_primer(no_ccdb_uha_dha_sgRNA_df,no_sgRNA_uha_dha_ccdb_df):

    #sgRNA质粒测序
    sgRNA_plasmid_sequencing_primer_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','plasmid','promoter_N20_terminator','promoter_N20_terminator_up','promoter_N20_terminator_down','seq_altered']]
    sgRNA_plasmid_sequencing_primer_df['plasmid_sequencing_region'] = sgRNA_plasmid_sequencing_primer_df['promoter_N20_terminator']
    sgRNA_plasmid_sequencing_primer_df['Region'] = sgRNA_plasmid_sequencing_primer_df['Name']+';'+sgRNA_plasmid_sequencing_primer_df['Region']
    sgRNA_plasmid_sequencing_primer_template = sgRNA_plasmid_sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
    sgRNA_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sgRNA_plasmid_sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')
    #ccdb质粒载体测序
    ccdb_plasmid_sequencing_primer_df = no_sgRNA_uha_dha_ccdb_df[['Name','Region','UHA','DHA','seq_altered','plasmid']]
    ccdb_plasmid_sequencing_primer_df['plasmid_sequencing_region'] = ccdb_plasmid_sequencing_primer_df['UHA'] + ccdb_plasmid_sequencing_primer_df['seq_altered'] + ccdb_plasmid_sequencing_primer_df['DHA']
    ccdb_plasmid_sequencing_primer_df['Region'] =  ccdb_plasmid_sequencing_primer_df['Name'] + ';' + ccdb_plasmid_sequencing_primer_df['Region']
    ccdb_plasmid_sequencing_primer_template = ccdb_plasmid_sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
    ccdb_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(ccdb_plasmid_sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')


    sgRNA_temp = no_ccdb_uha_dha_sgRNA_df[['Name','Region','promoter_N20_terminator']]
    sgRNA_temp['Region'] = sgRNA_temp['Name'] + ';' + sgRNA_temp['Region']
    ccdb_temp = no_sgRNA_uha_dha_ccdb_df[['Name','Region','UHA','DHA','seq_altered']]
    ccdb_temp['Region'] = ccdb_temp['Name'] + ';' + ccdb_temp['Region']
    sgRNA_plasmid_sequencing_primer_template = pd.merge(sgRNA_plasmid_sequencing_primer_template, sgRNA_temp, on='Region', how='inner')
    ccdb_plasmid_sequencing_primer_template = pd.merge(ccdb_plasmid_sequencing_primer_template, ccdb_temp,  on='Region', how='inner')

    return sgRNA_plasmid_sequencing_primer_df,sgRNA_plasmid_sequencing_primer_template, ccdb_plasmid_sequencing_primer_df,ccdb_plasmid_sequencing_primer_template

#genome sequencing primer
def genome_sequencing_design_primer(info_input_df, uha_dha_df):

    #编辑基因组测序
    info_input_df1 = info_input_df[['name','region','seq_uha_max_whole','seq_dha_max_whole','uha_upstream','dha_downstream']].rename(columns={'name':'Name','region':'Region'})
    UHA_DHA_df = pd.merge(info_input_df1,uha_dha_df)
    UHA_DHA_df['sequencing_template'] = UHA_DHA_df['uha_upstream'] + UHA_DHA_df['seq_uha_max_whole']+UHA_DHA_df['seq_altered']+UHA_DHA_df['seq_dha_max_whole'] + UHA_DHA_df['dha_downstream'] 
    UHA_DHA_df['sequencing_region'] = UHA_DHA_df['UHA']+UHA_DHA_df['seq_altered']+UHA_DHA_df['DHA']
    UHA_DHA_df['Region'] =  UHA_DHA_df['Name']+';'+UHA_DHA_df['Region']
    genome_sequencing_primer_df = p_d_seq.create_sequencing_primer(UHA_DHA_df,sr,'sequencing_template','sequencing_region',seq_type='genome_seq')
    #测序模板
    genome_sequencing_template = UHA_DHA_df[['Region','sequencing_template','UHA','seq_altered','DHA']]

    return genome_sequencing_primer_df, genome_sequencing_template
  

#     
def execute_one_plasmid_system(plasmid_primer_desgin_type,
                                sgRNA_region_seq_json,
                                sgRNA_primer_json,
                                gb_path,
                                info_df,
                                info_input_df,
                                uha_dha_df,
                                uha_dha_sgRNA_df,
                                uha_dha_info_primer_df,
                                uha_dha_primer_df,
                                enzyme_df,  
                                enzyme_name,
                                output
                                ):
    gb = SeqIO.read(gb_path, "genbank")
    gb_name = gb.name
    
    #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物  
    if plasmid_primer_desgin_type == 2:    
        #无需用户指定范围
        uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_p_df,type_kind = one_plasmid_system_pcr_design_primer(                                                                                                                            
                                                                                                                                        gb_path,
                                                                                                                                        info_df,
                                                                                                                                        uha_dha_sgRNA_df,
                                                                                                                                        uha_dha_info_primer_df,
                                                                                                                                        uha_dha_primer_df,
                                                                                                                                        enzyme_df,
                                                                                                                                        enzyme_name,
                                                                                                                                        plasmid_primer_desgin_type,
                                                                                                                                        sgRNA_region_seq_json
                                                                                                                                        )
        #修改 n20down_primer_p_df
        temp = n20down_primer_p_df.loc[:0]
        temp.loc[0,'Region'] = 1
        n20down_primer_p_df = temp

    elif plasmid_primer_desgin_type == 1:
        #用户指定范围
        uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_p_df,type_kind = one_plasmid_system_pcr_design_primer(                                                                                                                            
                                                                                                                                        gb_path,
                                                                                                                                        info_df,
                                                                                                                                        uha_dha_sgRNA_df,
                                                                                                                                        uha_dha_info_primer_df,
                                                                                                                                        uha_dha_primer_df,
                                                                                                                                        enzyme_df,
                                                                                                                                        enzyme_name,
                                                                                                                                        plasmid_primer_desgin_type,
                                                                                                                                        sgRNA_region_seq_json
                                                                                                                                        )
    elif plasmid_primer_desgin_type == 3:
         #用户指定额外引物  
        uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_p_df,type_kind = one_plasmid_system_pcr_design_primer(                                                                                                                            
                                                                                                                                        gb_path,
                                                                                                                                        info_df,
                                                                                                                                        uha_dha_sgRNA_df,
                                                                                                                                        uha_dha_info_primer_df,
                                                                                                                                        uha_dha_primer_df,
                                                                                                                                        enzyme_df,
                                                                                                                                        enzyme_name,
                                                                                                                                        plasmid_primer_desgin_type,
                                                                                                                                        sgRNA_primer_json
                                                                                                                                        )

    #设计质粒测序引物
    plasmid_sequencing_primer_df, sequencing_primer_template = one_plasmid_system_sequencing_design_primer(type_kind,uha_dha_sgRNA_df)

 
    #设计基因组测序引物
    genome_sequencing_primer_df, genome_sequencing_template = genome_sequencing_design_primer(info_input_df, uha_dha_df)

    #标准化，重命名   
    df_common_list = su.rename_common_primer_df(n20up_primer_p_df, n20down_primer_p_df, seq_altered_p_df)
    n20up_primer_p_df, n20down_primer_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1], df_common_list[2]
    uha_primer_df, dha_primer_df = su.rename_u_d_primer_df(uha_primer_df, dha_primer_df)
    df_sequencing_list = su.rename_sequencing_primer_df(plasmid_sequencing_primer_df, genome_sequencing_primer_df)
    plasmid_sequencing_primer_df, genome_sequencing_primer_df = df_sequencing_list[0], df_sequencing_list[1]

    #----------------------------生成gb文件用于可视化展示-------------------------------------------------------------------
    plasmid_primer_featrue_df = su.create_plasmid_primer_featrue_df(sequencing_primer_template,
                                                                 uha_primer_df,
                                                                 seq_altered_p_df,
                                                                 dha_primer_df,   
                                                                 n20up_primer_p_df)
    plasmid_primer_featrue_df = plasmid_primer_featrue_df.fillna('')
    joint_len, cut_seq_len = su.get_joint_by_enzyme(enzyme_df,enzyme_name)

    #为每个编辑区域创建gb文件
    pcr_gb_output = os.path.join(output,'one_plasmid_system_pcr_gb/')
    if not exists(pcr_gb_output):
        os.makedirs(pcr_gb_output)
    pcr_tsv_df = p_d_seq.create_gb_for_region(plasmid_primer_featrue_df, n20down_primer_p_df, joint_len, cut_seq_len, pcr_gb_output, type='sgRNA_ccdb')


    #质粒测序引物模板  
    plasmid_sequencing_template = sequencing_primer_template
    plasmid_sequencing_template.rename(columns={'Region':'ID','plasmid':'PLASMID','seq_altered':"SEQ_ALTERED",'promoter_N20_terminator':'PROMOTER_N20_TERMINATOR'},inplace=True)
    plasmid_seq_gb_output = os.path.join(output,'one_plasmid_system_plasmid_sequencing_gb/')
    if not exists(plasmid_seq_gb_output):
        os.makedirs(plasmid_seq_gb_output)
    plasmid_seq_tsv_df = p_d_seq.create_gb_for_sequencing_region(plasmid_sequencing_template, plasmid_sequencing_primer_df, plasmid_seq_gb_output, type='plasmid_sequencing')

    #基因组测序引物模板
    genome_sequencing_template.rename(columns={'Region':'ID','sequencing_template':'PLASMID','seq_altered':"SEQ_ALTERED"},inplace=True)
    genome_seq_gb_output = os.path.join(output,'one_plasmid_system_genome_sequencing_gb/')
    if not exists(genome_seq_gb_output):
        os.makedirs(genome_seq_gb_output)

    genome_seq_tsv_df = p_d_seq.create_gb_for_sequencing_region(genome_sequencing_template, genome_sequencing_primer_df, genome_seq_gb_output, type='genome_sequencing')

    #合并三个df 
    tsv_df = pd.merge(pcr_tsv_df,plasmid_seq_tsv_df,on='name',how='inner').merge(genome_seq_tsv_df,on='name',how='inner')
    tsv_df.to_csv(os.path.join(output,'one_plasmid_system_gb_visualization.tsv'), index=False, sep='\t') 
    #-----------------------------------------------------------------------------------------------------------------------

    #输出引物  
    xlsx_file = os.path.join(
        output,
        'one_plasmid_design_result.xlsx'
    )
    with pd.ExcelWriter(xlsx_file) as writer:  
        uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
        dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
        n20up_primer_p_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')  
        n20down_primer_p_df.to_excel(writer,sheet_name = 'Primer_plasmid',index_label='No.')
        seq_altered_p_df.to_excel(writer,sheet_name = 'Seq_altered',index_label='No.')
        plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
        genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')

    return xlsx_file                                                     

def execute_two_plasmid_system(
                                 method,
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
                                 sgRNA_region_seq_json,
                                 ccdb_region_seq_json,
                                 ccdb_label,
                                 promoter_terminator_label,
                                 n_20_label,
                                 output
                                 ):                                     

    no_ccdb_uha_dha_sgRNA_df, sgRNA_plasmid_backbone, promoter_seq, terminator_seq, sgRNA_promoter_terminator = p_d_seq.create_new_plasmid(no_ccdb_plasmid, uha_dha_sgRNA_df.copy(), ccdb_label, promoter_terminator_label, n_20_label)
    no_sgRNA_uha_dha_ccdb_df, ccdB_plasmid_backbone, ccdB_promoter_terminator_up_seq = p_d_seq.create_new_plasmid(no_sgRNA_plasmid, uha_dha_sgRNA_df.copy(), ccdb_label, promoter_terminator_label, n_20_label)
    #酶切退火方式
    

    #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物, 
    if  plasmid_primer_desgin_type == 1 and method == 'PCR':
        #对质粒进行分割
        #初步分割
        sgRNA_plasmid_seq, sgRNA_plasmid_region_seq = p_d_seq.plasmid_region_division_by_labels(no_ccdb_plasmid, ccdb_label, promoter_terminator_label, n_20_label)
        ccdb_plasmid_seq, ccdb_plasmid_region_seq = p_d_seq.plasmid_region_division_by_labels(no_sgRNA_plasmid, ccdb_label, promoter_terminator_label, n_20_label)

        #第一段距离：sgRNA终止子距离第一段区域的最小、最大距离
        #第二段距离：是第一段区域距离第二段区域的最小、最大距离
        #。。。。。。。。。
        #最后一段距离:最后一个区域距离sgRNA启动子距离


        #用户确定区域，程序给你设计模板，保证引物的数量和质量
        #将区域的序列转换成质粒上的坐标
        
        sgRNA_region_json = su.convert_seq_cor(no_ccdb_plasmid, sgRNA_region_seq_json)
        ccdb_region_json = su.convert_seq_cor(no_sgRNA_plasmid, ccdb_region_seq_json)

        #区域转换距离
        # sgRNA_region_json={
        #     "region1":"371,570",
        #     "region2":"3572,3770"
        # },
        # ccdb_region_json={
        #     "region1":"8364,8563",
        #     "region2":"376,575"
        # }
        # distance_dict = {
        #     'distance1':(3000,3200),
        #     'distance2':(3000,3200),
        # }
        #根据区域设计引物
        sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_region(
                                                                                    no_ccdb_uha_dha_sgRNA_df,
                                                                                    ccdB_plasmid_backbone,
                                                                                    no_sgRNA_plasmid,
                                                                                    no_ccdb_plasmid,
                                                                                    uha_dha_sgRNA_df,
                                                                                    enzyme_df,
                                                                                    enzyme_name,
                                                                                    sgRNA_plasmid_seq,
                                                                                    sgRNA_plasmid_region_seq,
                                                                                    ccdb_plasmid_seq,
                                                                                    ccdb_plasmid_region_seq,
                                                                                    sgRNA_region_json,
                                                                                    ccdb_region_json)

    elif    plasmid_primer_desgin_type == 2 and method == 'PCR':
            sgRNA_plasmid_p_df, ccdb_plasmid_p_df, sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df = two_plasmid_system_design_by_no_user(no_ccdb_uha_dha_sgRNA_df, ccdB_plasmid_backbone,enzyme_df,enzyme_name)
    elif    plasmid_primer_desgin_type == 3 and method == 'PCR':
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
        sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer(
                                                                                                    no_ccdb_uha_dha_sgRNA_df,
                                                                                                    uha_dha_sgRNA_df,
                                                                                                    sgRNA_primer_json, 
                                                                                                    ccdb_primer_json,
                                                                                                    sgRNA_plasmid_backbone,
                                                                                                    sgRNA_promoter_terminator,
                                                                                                    ccdB_plasmid_backbone,
                                                                                                    no_ccdb_plasmid,
                                                                                                    no_sgRNA_plasmid,
                                                                                                    enzyme_df,
                                                                                                    enzyme_name,
                                                                                                    n_20_label,
                                                                                                    promoter_terminator_label,
                                                                                                    ccdb_label)
    elif    plasmid_primer_desgin_type == 1 and method == 'OLIGO':
            #设计sgRNA质粒
            enzymeCutSeq_and_N20_df = two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df, promoter_seq, enzyme_df, enzyme_name)

            #设计ccdb质粒
            ccdb_plasmid_seq, ccdb_plasmid_region_seq = p_d_seq.plasmid_region_division_by_labels(no_sgRNA_plasmid, ccdb_label, promoter_terminator_label, n_20_label)
            ccdb_region_json = su.convert_seq_cor(no_sgRNA_plasmid, ccdb_region_seq_json)

            ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_ccdb_region(uha_dha_sgRNA_df,no_sgRNA_plasmid,ccdb_plasmid_seq,ccdb_plasmid_region_seq,ccdb_region_json, enzyme_df, enzyme_name)

    elif    plasmid_primer_desgin_type == 2 and method == 'OLIGO':
            #设计sgRNA质粒
            enzymeCutSeq_and_N20_df = two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df, promoter_seq, enzyme_df, enzyme_name)
            #设计ccdb质粒
            ccdb_plasmid_primer_df = two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name)

    elif    plasmid_primer_desgin_type == 3 and method == 'OLIGO':
            #设计sgRNA质粒
            enzymeCutSeq_and_N20_df = two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df, promoter_seq, enzyme_df, enzyme_name)
            #设计ccdb质粒
            ccdb_primer_position_json, ccdb_failture_primer = p_d_seq.check_locate_primer(ccdB_plasmid_backbone, ccdb_primer_json)

            if ccdb_primer_position_json == {} and ccdb_primer_json != {}:
                    return ccdb_failture_primer
            
            elif    ccdb_primer_position_json != {} and ccdb_primer_json != {}:
                    #设计ccdb primer
                    ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone,
                                                                                ccdb_primer_json,
                                                                                ccdb_primer_position_json,
                                                                                no_sgRNA_plasmid,
                                                                                enzyme_df,
                                                                                enzyme_name,
                                                                                n_20_label,
                                                                                ccdb_label,
                                                                                promoter_terminator_label)
            




    #设计seq_altered_primer，seq_altered > 120
    seq_altered_primer_template =  info_df[info_df.seq_altered.apply(lambda x:len(x)>120)][['Name','Region','seq_altered']]
    seq_altered_primer_template['Region'] = seq_altered_primer_template['Name'] +';'+ seq_altered_primer_template['Region']
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

    if method == 'PCR':
        sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df
    ccdb_plasmid_p_df = ccdb_plasmid_primer_df


    # if 'index' in sgRNA_plasmid_primer_df.columns and 'index' in ccdb_plasmid_primer_df.columns:
       
    #     sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df
    #     ccdb_plasmid_p_df = ccdb_plasmid_primer_df



    
    # else:
    #     if method =='PCR':
    #         sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
    #         ccdb_plasmid_p_df = ccdb_plasmid_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

    #设计质粒测序引物
    sgRNA_plasmid_sequencing_primer_df,sgRNA_plasmid_sequencing_primer_template, ccdb_plasmid_sequencing_primer_df,ccdb_plasmid_sequencing_primer_template = two_plasmid_system_sequencing_design_primer(no_ccdb_uha_dha_sgRNA_df,no_sgRNA_uha_dha_ccdb_df)

    #设计基因组测序引物
    genome_sequencing_primer_df,genome_sequencing_template = genome_sequencing_design_primer(info_input_df, uha_dha_df)

    #标准化，重命名
    if method == 'PCR':
        df_common_list = su.rename_common_primer_df(sgRNA_plasmid_p_df,ccdb_plasmid_p_df,seq_altered_p_df)
        sgRNA_plasmid_p_df, ccdb_plasmid_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1], df_common_list[2]
    else:
        df_common_list = su.rename_common_primer_df(ccdb_plasmid_p_df,seq_altered_p_df)
        ccdb_plasmid_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1]

    uha_primer_df, dha_primer_df = su.rename_u_d_primer_df(uha_primer_df,dha_primer_df)
    df_sequencing_list = su.rename_sequencing_primer_df(sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df,genome_sequencing_primer_df)
    sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df, genome_sequencing_primer_df = df_sequencing_list[0], df_sequencing_list[1], df_sequencing_list[2]

    
    
    #------------------------------------------生成gb文件用于引物的可视化展示-------------------------------------------------------------
    #--------------------------PCR----------------------生成sgRNA_gb文件-------------------------------------------------------------------
    sgRNA_plasmid_sequencing_primer_template = sgRNA_plasmid_sequencing_primer_template[['Region','plasmid','promoter_N20_terminator']].rename(columns={'Region':'ID','plasmid':"PLASMID","promoter_N20_terminator":"PROMOTER_N20_TERMINATOR"})
    if method == 'PCR':
        sgRNA_plasmid_primer = sgRNA_plasmid_p_df[['ID', 'PRIMER_LEFT_WHOLE_SEQUENCE', 'PRIMER_RIGHT_WHOLE_SEQUENCE']]
    else:
        sgRNA_plasmid_primer = enzymeCutSeq_and_N20_df[['ID','Target sequence']]
    joint_len, cut_seq_len = su.get_joint_by_enzyme(enzyme_df,enzyme_name)
    
    #为每个编辑区域创建gb文件
    gb_output = os.path.join(output,'two_plasmid_system_pcr_gb/')
    if not exists(gb_output):
        os.makedirs(gb_output)   
    if method == 'PCR':
        type = 'sgRNA'
        sgRNA_pcr_tsv_df = p_d_seq.create_gb_for_region(sgRNA_plasmid_sequencing_primer_template, sgRNA_plasmid_primer, joint_len, cut_seq_len, gb_output,type)
    elif method == 'OLIGO':
        type = 'enzyme_cut'
        sgRNA_pcr_tsv_df = p_d_seq.create_gb_for_region(sgRNA_plasmid_sequencing_primer_template, sgRNA_plasmid_primer, joint_len, cut_seq_len, gb_output,type)


    #---------------------------PCR---------------------生成ccdb_gb文件---------------------------------------------------------------------
    # plasmid_primer_featrue_df = ccdb_plasmid_sequencing_primer_template[['Region','plasmid']].rename(columns={'Region':'ID','plasmid':"PLASMID"})
    # [['Name','Region','UHA','DHA','seq_altered']]
    plasmid_primer_featrue_df = su.create_plasmid_primer_featrue_df(ccdb_plasmid_sequencing_primer_template,
                                                                 uha_primer_df,
                                                                 seq_altered_p_df,
                                                                 dha_primer_df)
    plasmid_primer_featrue_df = plasmid_primer_featrue_df.fillna('')
    joint_len, cut_seq_len = su.get_joint_by_enzyme(enzyme_df,enzyme_name)

    #为每个编辑区域创建gb文件
    gb_output = os.path.join(output ,'two_plasmid_system_pcr_gb/')
    if not exists(gb_output):
        os.makedirs(gb_output)   
    ccdb_pcr_tsv_df = p_d_seq.create_gb_for_region(plasmid_primer_featrue_df, ccdb_plasmid_p_df,joint_len, cut_seq_len, gb_output,type='ccdb')




















    #---------------------plasmid-------SEQUENCING------------------------------------------------------------------------------------------------------------
    #为每个编辑区域创建gb文件
    gb_output = os.path.join(output,'two_plasmid_system_plasmid_sequencing_gb/')
    if not exists(gb_output):
        os.makedirs(gb_output)
    sgRNA_plasmid_primer_featrue_df = sgRNA_plasmid_sequencing_primer_template
    sgRNA_plasmid_sequencing_tsv_df = p_d_seq.create_gb_for_sequencing_region(sgRNA_plasmid_primer_featrue_df, sgRNA_plasmid_sequencing_primer_df, gb_output, type='sgRNA_plasmid_sequencing')


    ccdb_plasmid_primer_featrue_df = ccdb_plasmid_sequencing_primer_template[['Region','plasmid','UHA','DHA']].rename(columns={'Region':'ID','plasmid':"PLASMID"})
    ccdb_plasmid_sequencing_tsv_df = p_d_seq.create_gb_for_sequencing_region(ccdb_plasmid_primer_featrue_df, ccdb_plasmid_sequencing_primer_df, gb_output, type='ccdb_plasmid_sequencing')

    #-------------------genome----------SEQUENCING------------------------------------------------------------------
    gb_output = os.path.join(output,'two_plasmid_system_genome_sequencing_gb/')
    if not exists(gb_output):
        os.makedirs(gb_output)
    genome_sequencing_template.rename(columns={'Region':'ID','sequencing_template':'PLASMID','seq_altered':"SEQ_ALTERED"},inplace=True)
    genome_sequencing_tsv_df = p_d_seq.create_gb_for_sequencing_region(genome_sequencing_template, genome_sequencing_primer_df, gb_output, type='genome_sequencing')


    #生成tsv
    # print(sgRNA_tsv_df,'\n',ccdb_tsv_df)  
    pcr_df = pd.merge(sgRNA_pcr_tsv_df, ccdb_pcr_tsv_df, on='name')  
    sequencing_df = pd.merge(sgRNA_plasmid_sequencing_tsv_df, ccdb_plasmid_sequencing_tsv_df,on='name')
    tsv_df = pcr_df.merge(sequencing_df,on='name').merge(genome_sequencing_tsv_df,on='name')
    tsv_df.to_csv(os.path.join(output,'two_plasmid_system_gb_visualization.tsv'), index=False, sep='\t')
    #-------------------------------------------------------------------------------------------------------------------------------------------


    #输出引物     
    xlsx_file = os.path.join(
        output,
        'two_plasmid_design_result.xlsx'
    )
    if method == 'PCR':
        with pd.ExcelWriter(xlsx_file) as writer:  
            uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            sgRNA_plasmid_p_df.to_excel(writer,sheet_name = 'sgRNA_plasmid_fragment',index_label='No.')  
            ccdb_plasmid_p_df.to_excel(writer,sheet_name = 'uha_dha_plasmid',index_label='No.')
            seq_altered_p_df.to_excel(writer,sheet_name = 'Seq_altered',index_label='No.')
            sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
            ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
            genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')
    else:
         with pd.ExcelWriter(xlsx_file) as writer:  
            uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            enzymeCutSeq_and_N20_df.to_excel(writer,sheet_name = 'sgRNA_plasmid_fragment',index_label='No.')  
            ccdb_plasmid_p_df.to_excel(writer,sheet_name = 'uha_dha_plasmid',index_label='No.')
            seq_altered_p_df.to_excel(writer,sheet_name = 'Seq_altered',index_label='No.')
            sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
            ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
            genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')

    return xlsx_file

  
def read_chopchopInput_add_uha_dha(genome_path,chopchop_input,uha_dha_params):
    max_left_arm_seq_length = uha_dha_params['max_left_arm_seq_length']
    max_right_arm_seq_length = uha_dha_params['max_right_arm_seq_length']

    info_input_df = su.del_Unnamed(pd.read_csv(chopchop_input))

    info_input_df.columns = [i.lower() for i in info_input_df.columns]


    def work(mun_id,geneid, mutation_pos_index):
        if mutation_pos_index - max_left_arm_seq_length < 0:
            error_message = "The length of upstream sequence of manipulation site of " + mun_id + " must be larger than sum of 'Max Length of UHA' and 'Max Length of UIS'."
            return error_message,error_message,error_message,error_message

        record = su.extract_seq_from_genome(genome_path,geneid)

        seq_uha_max_whole = str(record[
                        mutation_pos_index - max_left_arm_seq_length : mutation_pos_index
                        ])
        seq_dha_max_whole = str(record[
                                mutation_pos_index : mutation_pos_index + max_right_arm_seq_length
                                ])


        uha_upstream = str(  
                        record[
                            mutation_pos_index - max_left_arm_seq_length - 100 : mutation_pos_index - max_left_arm_seq_length
                        ]
                    )
        dha_downstream=str(
                        record[
                            mutation_pos_index + max_right_arm_seq_length  : mutation_pos_index + max_right_arm_seq_length  + 100
                        ]
                    )
        return  uha_upstream, dha_downstream, seq_uha_max_whole, seq_dha_max_whole

    info_df = su.lambda2cols(info_input_df,lambdaf=work,in_coln=['name','geneid','mutation_pos_index'],to_colns=['uha_upstream','dha_downstream','seq_uha_max_whole','seq_dha_max_whole'])

    return info_df  


def check_quality_control(gb_path,seq_json):
    seq_json = su.check_seq_in_gb(gb_path, seq_json)  
    failture_seq_json = {}
    for k,v in seq_json.items():
        if 'The sequence' in v:
            failture_seq_json.update({k:v})

    return failture_seq_json,seq_json

   
def check_plasmid(gb_path, ccdb_label='', promoter_terminator_label='', n_20_label=''):

    gb = SeqIO.read(gb_path, "genbank")
    gb_seq = str(gb.seq)
    #get coordinate
    ccdb_coordinate = su.get_feature_coordinate(ccdb_label,gb_path) 
    #N20
    n20_coordinate = su.get_feature_coordinate(n_20_label,gb_path)

    if ccdb_coordinate != (-1,-1) and n20_coordinate !=(-1,-1):
        return 'one_plasmid_file_path'
    elif ccdb_coordinate==(-1,-1) and n20_coordinate !=(-1,-1):
        return 'no_ccdb_plasmid'
    elif ccdb_coordinate != (-1,-1) and n20_coordinate ==(-1,-1):
        return 'no_sgRNA_plasmid'
    else:
        return 'error'


def check_enzyme(enzyme,enzyme_df):
    name = enzyme['enzyme_name']
    gap_seq = enzyme['gap_sequence']
    protection_seq = enzyme['protection_sequence']

    if name in list(enzyme_df['name']) and len(enzyme_df[enzyme_df['name']==name]) == 1:
        enzyme_df = enzyme_df[enzyme_df['name']==name]
        enzyme_df['protective_base'] = protection_seq
        enzyme_df['gap_len'] = len(protection_seq)
        enzyme_df['gap_seq'] = gap_seq
        return enzyme_df
    else:
        return 'There is a problem with the enzyme you provided'


def main(data):
    
    chopchop_input = data['chopchop_input']

    #uha_dha参数
    uha_dha_params = data['uha_dha_config']
    
    # enzyme_path = parent_base_path +'/'+ data['enzyme_path']
    sgRNA_result_path = data['sgRNA_result_path']
    plasmid_file_1 = data['one_plasmid_file_path']
    plasmid_file_2 = data['no_ccdb_plasmid']
    plasmid_file_3 = data['no_sgRNA_plasmid']
    
    #plasmid label
    ccdb_label = data['plasmid_label']['ccdb_label']
    promoter_terminator_label = data['plasmid_label']['promoter_terminator_label']
    n_20_label = data['plasmid_label']['n_20_label']   

    one_plasmid_file_path=''
    no_ccdb_plasmid=''
    no_sgRNA_plasmid=''
    
    #检查质粒,自动判断质粒类型
    if plasmid_file_1 != '':
        plasmid_type = check_plasmid(plasmid_file_1, ccdb_label, promoter_terminator_label, n_20_label)
        if plasmid_type == 'one_plasmid_file_path':
            one_plasmid_file_path = plasmid_file_1
        elif plasmid_type == 'no_ccdb_plasmid':
            no_ccdb_plasmid = plasmid_file_1
        elif plasmid_type == 'no_sgRNA_plasmid':
            no_sgRNA_plasmid = plasmid_file_1
        elif plasmid_type == 'error':
            return  'There is a problem with the plasmid you uploaded'
    
    if plasmid_file_2 != '':
        plasmid_type = check_plasmid(plasmid_file_2, ccdb_label, promoter_terminator_label, n_20_label)
        if plasmid_type == 'one_plasmid_file_path':
            one_plasmid_file_path = plasmid_file_2
        elif plasmid_type == 'no_ccdb_plasmid':
            no_ccdb_plasmid = plasmid_file_2
        elif plasmid_type == 'no_sgRNA_plasmid':
            no_sgRNA_plasmid = plasmid_file_2
        elif plasmid_type == 'error':
            return  'There is a problem with the plasmid you uploaded'

    if plasmid_file_3 != '':
        plasmid_type = check_plasmid(plasmid_file_3, ccdb_label, promoter_terminator_label, n_20_label)    
        if plasmid_type == 'one_plasmid_file_path':
            one_plasmid_file_path = plasmid_file_3
        elif plasmid_type == 'no_ccdb_plasmid':
            no_ccdb_plasmid = plasmid_file_3
        elif plasmid_type == 'no_sgRNA_plasmid':
            no_sgRNA_plasmid = plasmid_file_3
        elif plasmid_type == 'error':
            return  'There is a problem with the plasmid you uploaded'


    #primer
    sgRNA_primer_json = data['sgRNA_primer_json']
    ccdb_primer_json = data['ccdb_primer_json']

    primer_json = data['primer_json']  
    region_seq_json = data['region_json']


    if primer_json != {} and one_plasmid_file_path !='':
        failture_seq_json, seq_json = check_quality_control(one_plasmid_file_path ,primer_json)
        if failture_seq_json == {}:
            primer_json = seq_json
        else:
            return failture_seq_json
    
    if region_seq_json != {} and one_plasmid_file_path !='':
        failture_seq_json, seq_json = check_quality_control(one_plasmid_file_path ,region_seq_json)
        if failture_seq_json == {}:
            region_seq_json = seq_json
        else:
            return failture_seq_json

    if sgRNA_primer_json !={} and no_ccdb_plasmid !='':
        failture_seq_json, seq_json = check_quality_control(no_ccdb_plasmid,sgRNA_primer_json)
        if failture_seq_json == {}:
            sgRNA_primer_json = seq_json
        else:
            return failture_seq_json

    if  ccdb_primer_json !={} and no_sgRNA_plasmid != '':
        failture_seq_json, seq_json = check_quality_control(no_sgRNA_plasmid,ccdb_primer_json)
        if failture_seq_json == {}:
            ccdb_primer_json = seq_json
        else:
            return failture_seq_json

    # sgRNA_region_json = data['sgRNA_region_json']
    # ccdb_region_json = data['ccdb_region_json']

    sgRNA_region_seq_json = data['sgRNA_region_json']
    if  sgRNA_region_seq_json !={} and no_ccdb_plasmid !='':
        failture_seq_json, seq_json = check_quality_control(no_ccdb_plasmid,sgRNA_region_seq_json)
        if failture_seq_json == {} :
            sgRNA_region_seq_json = seq_json
        else:
            return failture_seq_json

    ccdb_region_seq_json = data['ccdb_region_json']
    if  sgRNA_region_seq_json !={} and no_sgRNA_plasmid !='':
        failture_seq_json, seq_json = check_quality_control(no_sgRNA_plasmid,ccdb_region_seq_json)
        if failture_seq_json == {}:
            ccdb_region_seq_json = seq_json  
        else:
            return failture_seq_json


    #genome
    genome_path = data['ref_genome']

    #配置引物参数
    # config.S_GLOBAL_ARGS = data['S_GLOBAL_ARGS']
    # config.Q_ARGS = data['Q_ARGS']
    config.UHA_ARGS = data['UHA_ARGS']
    config.SEQ_ALTERED_ARGS = data['SEQ_ALTERED_ARGS']
    config.DHA_ARGS = data['DHA_ARGS']

    if data.get('UP_SGRNA_ARGS') == {} or data.get('DOWN_SGRNA_ARGS') == {}:
        method = 'OLIGO'
        print('执行方法：', method)
    else:
        method = 'PCR'
        print('执行方法：', method)
        config.UP_SGRNA_ARGS = data['UP_SGRNA_ARGS']
        config.DOWN_SGRNA_ARGS = data['DOWN_SGRNA_ARGS']

    config.PLASMID_Q_ARGS = data['PLASMID_Q_ARGS']
    config.GENOME_Q_ARGS = data['GENOME_Q_ARGS']   
    
    #配置输出参数
    output = data['edit_sequence_design_workdir']
    if not os.path.exists(output):
        os.makedirs(output)

    scene = data['scene']  

    
    #0 参数的质量控制
    # 1.read 编辑序列信息,给chopchop输入加uha、dha信息
    info_input_df = read_chopchopInput_add_uha_dha(genome_path, chopchop_input, uha_dha_params)


    # 2.#读取用户填的酶
    # enzyme={
    #     "enzyme_name":"BsaI",  
    #     "protective_base":"CCA",
    #     "recognition_seq":"GGTCTC",
    #     "cut_seq_len":4,
    #     "gap_len":1    
    # }
    base_path = os.path.abspath(os.path.dirname(__file__)) + '/'
    enzyme_path =  base_path + '/input/enzyme.csv'      
    enzyme_df = su.del_Unnamed(pd.read_csv(enzyme_path))
    enzyme = data['enzyme']
    enzyme_name = enzyme['enzyme_name']     
    #检查酶的相关信息
    return_value = check_enzyme(enzyme, enzyme_df)
    if type(return_value) == str:
        return return_value
    else:
        enzyme_df = return_value




    # 3.提取用户选择的sgRNA     
    if scene == 'only_primer': 
        sgRNA = info_input_df[['name','crrna','region']].rename(columns={'crrna':"Target sequence","name":"Name","region":"Region"})
        sgRNA['Rev Target sequence'] = sgRNA['Target sequence'].apply(lambda x: su.revComp(x)) 

    elif scene == 'both_sgRNA_primer':
        selected_sgRNA_result = data['sgRNA_result']
        sgRNA = p_d_seq.extract_sgRNA_from_chopchop(sgRNA_result_path, selected_sgRNA_result)
    

    # 4.设计源生同源臂引物
    uha_dha_primer_df = extract_uha_dha_primer(info_input_df, sgRNA)   

    # 5.提取同源臂
    uha_dha_info_primer_df, uha_dha_df, uha_dha_sgRNA_df, info_df = extract_uha_dha(info_input_df, uha_dha_primer_df, sgRNA)
  
    #判断质粒的执行类型    
    if  one_plasmid_file_path != '' and no_ccdb_plasmid == '' and no_sgRNA_plasmid == '':
        plasmid_system_type = 1
    elif  one_plasmid_file_path =='' and no_ccdb_plasmid != '' and no_sgRNA_plasmid != '':
        plasmid_system_type = 2
    elif  one_plasmid_file_path !='' and no_ccdb_plasmid !='' and no_sgRNA_plasmid !='':
        plasmid_system_type = 0
    else:
        return '你选择这的双质粒系统，质粒没有上传完整!'
    
    print('--1.执行单质粒系统,--2.执行双质粒系统,---0.执行单、双质粒系统都执行------现在正在执行的情况：',plasmid_system_type) 

    if plasmid_system_type == 1:

        #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物
        if region_seq_json == {} and  primer_json == {}:
            plasmid_primer_desgin_type = 2
        elif region_seq_json != {} and primer_json == {}:
            plasmid_primer_desgin_type = 1
        elif primer_json != {} and region_seq_json == {}:
            plasmid_primer_desgin_type = 3
        else:
            plasmid_primer_desgin_type = 1

        # 6.执行单质粒系统
        one_plasmid_output_path = execute_one_plasmid_system(   
                                                                plasmid_primer_desgin_type,
                                                                region_seq_json,
                                                                primer_json,
                                                                one_plasmid_file_path,
                                                                info_df,
                                                                info_input_df,
                                                                uha_dha_df,
                                                                uha_dha_sgRNA_df,
                                                                uha_dha_info_primer_df,
                                                                uha_dha_primer_df,
                                                                enzyme_df,
                                                                enzyme_name,
                                                                output
                                                                )
        return one_plasmid_output_path
    
    elif plasmid_system_type == 2:
        # 7.执行双质粒系统
        #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物
       
        if sgRNA_primer_json !={} or ccdb_primer_json !={}:
            plasmid_primer_desgin_type = 3
        elif sgRNA_region_seq_json !={} or ccdb_region_seq_json !={}:
            plasmid_primer_desgin_type = 1
        elif sgRNA_primer_json == {} and ccdb_primer_json == {} and sgRNA_region_seq_json == {} and ccdb_region_seq_json == {}:
            plasmid_primer_desgin_type = 2

        print('--1---用户指定范围,--2----无需用户指定范围,--3----用户指定额外引物-----4.酶切酶连方式-------现在执行的情况：',plasmid_primer_desgin_type)    

        two_plasmid_output_path = execute_two_plasmid_system(
                                    method,
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
                                    sgRNA_region_seq_json,
                                    ccdb_region_seq_json,
                                    ccdb_label,  
                                    promoter_terminator_label,
                                    n_20_label,
                                    output
                                    )
        return two_plasmid_output_path

    elif plasmid_system_type == 0:

         
        #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物
        if region_seq_json == {} and  primer_json == {}:
            plasmid_primer_desgin_type = 2
        elif region_seq_json != {} and primer_json == {}:
            plasmid_primer_desgin_type = 1
        elif primer_json != {} and region_seq_json == {}:  
            plasmid_primer_desgin_type = 3
        elif primer_json != {} and region_seq_json != {}:
            plasmid_primer_desgin_type = 1
        print('--1.用户指定范围---2.无需用户指定范围---3.用户指定额外引物------现在正在执行的情况：',plasmid_primer_desgin_type)

        one_plasmid_output_path = execute_one_plasmid_system(   
                                                                plasmid_primer_desgin_type,
                                                                region_seq_json,
                                                                primer_json,
                                                                one_plasmid_file_path,
                                                                info_df,
                                                                info_input_df,
                                                                uha_dha_df,
                                                                uha_dha_sgRNA_df,
                                                                uha_dha_info_primer_df,
                                                                uha_dha_primer_df,
                                                                enzyme_df,
                                                                enzyme_name,
                                                                output
                                                                )

        
         # 7.执行双质粒系统
        #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物
       
        if sgRNA_primer_json !={} or ccdb_primer_json !={}:
            plasmid_primer_desgin_type = 3
        elif sgRNA_region_seq_json !={} or ccdb_region_seq_json !={}:
            plasmid_primer_desgin_type = 1
        elif sgRNA_primer_json == {} and ccdb_primer_json == {} and sgRNA_region_seq_json == {} and ccdb_region_seq_json == {}:
            plasmid_primer_desgin_type = 2

        print('--------------------------------------现在执行的情况：',plasmid_primer_desgin_type)    

        two_plasmid_output_path = execute_two_plasmid_system(
                                    method,
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
                                    sgRNA_region_seq_json,
                                    ccdb_region_seq_json,
                                    ccdb_label,  
                                    promoter_terminator_label,
                                    n_20_label,
                                    output
                                    )
        
    return one_plasmid_output_path, two_plasmid_output_path

     
if __name__ == '__main__':

    #read json
    # import argparse
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--input', '-i', help='input params file', required=True) 
    # args = parser.parse_args()
    # input_path =  args.input
    # with open(input_path, "r") as f:
    #     data = json.load(f)     

    # main(data)         
 
    data1 = {     
        "chopchop_input": "/home/yanghe/tmp/data_preprocessing/output/info_input.csv",   
        "sgRNA_result_path": "/home/yanghe/tmp/chopchop/output/sgRNA.csv",
        "edit_sequence_design_workdir":"/home/yanghe/tmp/edit_sequence_design/output/",
        "ref_genome":"/home/yanghe/program/data_preprocessing/input/GCA_000011325.1_ASM1132v1_genomic.fna",
        "one_plasmid_file_path":"/home/yanghe/program/edit_sequence_design/input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
        "no_ccdb_plasmid":"/home/yanghe/program/edit_sequence_design/input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "no_sgRNA_plasmid":"/home/yanghe/program/edit_sequence_design/input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "scene":"both_sgRNA_primer",
        "uha_dha_config": {
            "max_right_arm_seq_length": 1050,  
            "max_left_arm_seq_length": 1050,   
            "min_left_arm_seq_length": 1000,   
            "min_right_arm_seq_length": 1000     
        },

        "plasmid_label":{
            "ccdb_label":"ccdB",  
            "promoter_terminator_label":"gRNA",
            "n_20_label":"N20"
        },

        "primer_json":{
        
        },
        "region_json":{
            
        },     

        "sgRNA_primer_json":{
           
        },
        "ccdb_primer_json":{
                
        },   
    
        "sgRNA_region_json":{
              
        },
        
        "ccdb_region_json":{  
            
        },   
        
        "enzyme":{
            "enzyme_name":"BsaI",
            "gap_sequence":"A",  
            "protection_sequence":"CCA"   
        },  
        
        "UHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "SEQ_ALTERED_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,  
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "DHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },

        "PLASMID_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "GENOME_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,     
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        'sgRNA_result':{
            "Cgl0006_1176_G_A_sub":"1",
            "Cgl2342_213_GCA_ins":"1",
            "Cgl1436_1113_CAA_del":"1",
            "Cgl1790_1647_TCC_sub":"1",
            "Cgl1386_327_18to15_sub":"1",
            "Cgl0591_-1_Ppgk_promoter_ins":"1",
            "Cgl0141_cds_del":"1",
            "153019_ecoil_ybeL_ins":"1",
            "Cgl0851_ecoli_pgi_sub":"1"
        }      
    }


    data2 = {     
        "chopchop_input": "./input/only_primer/info_input.csv",   
        "sgRNA_result_path": "",
        "edit_sequence_design_workdir":"./output/only_primer",
        "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
        "one_plasmid_file_path":"./input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
        "no_ccdb_plasmid":"./input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "no_sgRNA_plasmid":"./input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "scene":"only_primer",
        "uha_dha_config": {
            "max_right_arm_seq_length": 1050,  
            "max_left_arm_seq_length": 1050,   
            "min_left_arm_seq_length": 1000,   
            "min_right_arm_seq_length": 1000     
        },

        "plasmid_label":{
            "ccdb_label":"ccdB",  
            "promoter_terminator_label":"gRNA",
            "n_20_label":"N20"
        },

        "primer_json":{
        
        },
        "region_json":{
            
        },     

        "sgRNA_primer_json":{
           
        },
        "ccdb_primer_json":{
                
        },   
    
        "sgRNA_region_json":{
              
        },
        
        "ccdb_region_json":{  
            
        },   
        
        "enzyme":{
            "enzyme_name":"BsaI",
            "gap_sequence":"A",  
            "protection_sequence":"CCA"   
        },  
        
        "UHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "SEQ_ALTERED_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,  
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "DHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },

        "PLASMID_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "GENOME_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,     
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "UP_SGRNA_ARGS": {
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 65,
            "PRIMER_MIN_GC": 30,
            "PRIMER_MAX_GC": 70
        },
        "DOWN_SGRNA_ARGS": {
           
        },
        'sgRNA_result':{
            "Cgl0006_1176_G_A_sub":"1",
            "Cgl2342_213_GCA_ins":"1",
            "Cgl1436_1113_CAA_del":"1",
            "Cgl1790_1647_TCC_sub":"1",
            "Cgl1386_327_18to15_sub":"1",
            "Cgl0591_-1_Ppgk_promoter_ins":"1",
            "Cgl0141_cds_del":"1",
            "153019_ecoil_ybeL_ins":"1",
            "Cgl0851_ecoli_pgi_sub":"1"
        }      
    }
    main(data2)     




        