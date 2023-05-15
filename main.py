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

    uha_dha_primer_df,failture_uha_dha_primer_df = p_d_seq.design_primer(primer_template_for_u_d_ha_df,'Name_Region','primer_template','u_d')
    uha_dha_primer_df = uha_dha_primer_df.rename(columns={'Region':'id'})
    uha_dha_primer_df = uha_dha_primer_df.join(uha_dha_primer_df.id.str.split(';',expand=True).rename(columns = {0:'Name',1:'Region',2:'Type'})).drop(columns='id')
    
    return uha_dha_primer_df, failture_uha_dha_primer_df

#uha_dha 
def extract_uha_dha(info_input_df, uha_dha_primer_df, sgRNA):
    #整合进突变信息
    info_df = info_input_df[['name','region','seq_altered','type','ref','strand']].rename(columns={'name':'Name','region':'Region'})
    info_df.seq_altered.fillna('',inplace=True)
    uha_dha_info_primer_df = pd.merge(info_df,uha_dha_primer_df,on=['Name','Region'], how='inner')

    #提取源生同源臂
    uha_dha_df = p_d_seq.create_uha_dha_df(uha_dha_primer_df) 
    #合并突变信息
    uha_dha_df = pd.merge(uha_dha_df,info_df,on=['Name','Region'])
    uha_dha_sgRNA_df = pd.merge(uha_dha_df,sgRNA,on=['Name','Region'],how='inner')

    return uha_dha_info_primer_df, uha_dha_df, uha_dha_sgRNA_df,info_df

def one_plasmid_system_design_by_user_primer(primer_json,primer_position_json,gb_path,plasmid_backbone,enzyme_df,enzyme_name):

    sgRNA_promoter_terminator_start = 0
    #对引物进行排序,确定所有引物
   
    state, primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         primer_json,
                                                         primer_position_json,
                                                         gb_path,
                                                         plasmid_backbone,
                                                         plasmid_type='one')
    plasmid_primer_df = pd.DataFrame()
    failture_plasmid_primer_df = pd.DataFrame()
    if state == 'success':    
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
    elif state == 'failture':
        failture_plasmid_primer_df = primer_dict_df       

    return plasmid_primer_df, failture_plasmid_primer_df

def one_plasmid_system_design_by_user_region(uha_dha_sgRNA_df,sgRNA_region_seq_json,one_plasmid_file_path,sgRNA_plasmid_backbone,enzyme_df,enzyme_name):
    
    plasmid_primer_df = pd.DataFrame()
    failture_plasmid_primer_df = pd.DataFrame()

    #序列转换成坐标
    region_json = su.convert_seq_cor(one_plasmid_file_path, sgRNA_region_seq_json,seq=sgRNA_plasmid_backbone)
    #坐标转转换成距离
    distance_dict = region_2_distance(len(sgRNA_plasmid_backbone), region_json, 0)   

    #根据区域设计引物
    state,primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(0,sgRNA_plasmid_backbone,distance_dict)
    if state == 'success':
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
       
    elif state == 'failture':
        failture_plasmid_primer_df = pd.DataFrame(primer_result_list)

    return plasmid_primer_df, failture_plasmid_primer_df

#one plasmid pcr primer    
def one_plasmid_system_pcr_design_primer(gb_path,
                                         info_df,
                                         uha_dha_sgRNA_df,
                                         uha_dha_info_primer_df,
                                         uha_dha_primer_df,
                                         enzyme_df,
                                         enzyme_name,
                                         plasmid_primer_desgin_type,
                                         seq_json,
                                         ccdb_label,
                                         promoter_terminator_label,
                                         n_20_label,
                                         promoter_label):
    
    #创建新的质粒      
    uha_dha_sgRNA_df,promoter_seq, promoter_terminator_up_promoter_seq, promoter_terminator_down_terminator_seq, type_kind  = p_d_seq.create_new_plasmid(gb_path, uha_dha_sgRNA_df.copy(), ccdb_label, promoter_terminator_label, n_20_label, promoter_label)
    uha_dha_sgRNA_df['promoter_seq'] = promoter_seq


    #设计sgRNA、ccdb质粒引物
    n20up_primer_template = uha_dha_sgRNA_df[['Name','Region','n20_up_template','Target sequence','Rev Target sequence']]
    n20up_primer_template['Region'] = n20up_primer_template['Name'] +';'+ n20up_primer_template['Region']
    n20up_primer_df, failture_n20up_primer_df = p_d_seq.design_primer(n20up_primer_template,'Region','n20_up_template','sgRNA')
    if len(n20up_primer_df)>0:
        n20up_primer_df = pd.merge(n20up_primer_template[['Region','Target sequence','Rev Target sequence']],n20up_primer_df,on=['Region'],how='inner')

        #ccdb、sgrna质粒引物加接头
        n20up_primer_df = p_d_seq.add_joint_sgRNA_primer(n20up_primer_df,enzyme_df, enzyme_name,'',stype='n20up_primer_joint')


     #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物  
    if plasmid_primer_desgin_type == 2: 
        n20down_primer_template = uha_dha_sgRNA_df[['Name','Region','n20_down_template','Target sequence','Rev Target sequence']]
        n20down_primer_template['Region'] = n20down_primer_template['Name'] +';'+ n20down_primer_template['Region']
        n20down_primer_df, failture_n20down_primer_df = p_d_seq.design_primer(n20down_primer_template,'Region','n20_down_template','sgRNA')
        n20down_primer_df = pd.merge(n20down_primer_template[['Region','Target sequence','Rev Target sequence']],n20down_primer_df,on=['Region'],how='inner')
        #加接头
        n20down_primer_df = p_d_seq.add_joint_sgRNA_primer(n20down_primer_df,enzyme_df,enzyme_name,'',stype='n20down_primer_joint')

    elif plasmid_primer_desgin_type == 1:
        plasmid_backbone = promoter_terminator_down_terminator_seq
        n20down_primer_df, failture_n20down_primer_df =  one_plasmid_system_design_by_user_region(uha_dha_sgRNA_df,seq_json,gb_path,plasmid_backbone,enzyme_df,enzyme_name)

    elif plasmid_primer_desgin_type == 3:
        plasmid_backbone = promoter_terminator_down_terminator_seq
        primer_position_json, sgRNA_failture_primer = p_d_seq.check_locate_primer(plasmid_backbone, seq_json)
        n20down_primer_df, failture_n20down_primer_df = one_plasmid_system_design_by_user_primer(seq_json, primer_position_json, gb_path, plasmid_backbone, enzyme_df, enzyme_name)


     #给同源臂引物加接头:uha取promoter_terminator_down_terminator_seq尾部反义4bp，dha取头promoter_terminator_up_promoter_seq正义4bp
    uha_dha_primer_df = p_d_seq.add_joint_sgRNA_primer(uha_dha_info_primer_df,enzyme_df,enzyme_name,promoter_terminator_down_terminator_seq, promoter_terminator_up_promoter_seq, stype = 'u_d_primer_joint')

    seq_altered_primer_template = info_df[info_df.seq_altered.apply(lambda x:len(x)>120)][['Name','Region','seq_altered']]

    if len(seq_altered_primer_template)>0:
        seq_altered_primer_template['Region'] = seq_altered_primer_template['Name'] +';'+ seq_altered_primer_template['Region']
        seq_altered_primer_df, failture_seq_altered_primer_df  = p_d_seq.design_primer(seq_altered_primer_template,'Region','seq_altered',stype='seq_altered')
        #seq_altered_primer加接头
        if len(seq_altered_primer_df) > 0:
            seq_altered_primer_df = p_d_seq.add_joint_sgRNA_primer(seq_altered_primer_df,enzyme_df,enzyme_name,'',stype='seq_altered_primer_joint')
        else:
            seq_altered_primer_df = pd.DataFrame()
    else:
        seq_altered_primer_df = pd.DataFrame()
        failture_seq_altered_primer_df = pd.DataFrame()

    #分别提取，修饰后的uha、dha引物
    in_col = ['Name','Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']
    ou_col = ['Name','Region',"u_primer_f_seq_(5'-3')","u_primer_r_seq_(5'-3')",'UHA','UHA_size',"d_primer_f_seq_(5'-3')","d_primer_r_seq_(5'-3')",'DHA','DHA_size']
    uha_primer_df, dha_primer_df = p_d_seq.create_uha_dha_primer_df(uha_dha_primer_df, in_col, ou_col)

    n20down_primer_p_df = pd.DataFrame()
    if len(n20down_primer_df)>0:
        n20down_primer_p_df = n20down_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]
    n20up_primer_p_df = pd.DataFrame() 
    if len(n20up_primer_df)>0:
        n20up_primer_p_df = n20up_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]

    if len(seq_altered_primer_df)>0:
        seq_altered_primer_df = seq_altered_primer_df[["Region","primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']]

    return uha_dha_sgRNA_df, uha_primer_df, dha_primer_df, n20down_primer_p_df, n20up_primer_p_df, seq_altered_primer_df,failture_seq_altered_primer_df, type_kind, failture_n20up_primer_df, failture_n20down_primer_df
   
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
        plasmid_sequencing_primer_df1, failture_plasmid_sequencing_primer_df1 = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')
        #promoter_N20_terminator测序
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['promoter_N20_terminator']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df2, failture_plasmid_sequencing_primer_df2 = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')
        
        if len(plasmid_sequencing_primer_df1) >0 and len(plasmid_sequencing_primer_df2)>0:
            plasmid_sequencing_primer_df = p_d_seq.merge_sequencing_result(plasmid_sequencing_primer_df1, plasmid_sequencing_primer_df2)
        elif len(plasmid_sequencing_primer_df1) > 0 and len(plasmid_sequencing_primer_df2)==0: 
            plasmid_sequencing_primer_df = plasmid_sequencing_primer_df1
        elif len(plasmid_sequencing_primer_df2) > 0 and len(plasmid_sequencing_primer_df1) == 0:
            plasmid_sequencing_primer_df = plasmid_sequencing_primer_df2
        
        failture_plasmid_sequencing_primer_df = failture_plasmid_sequencing_primer_df1.append(failture_plasmid_sequencing_primer_df2) 

    elif type_kind == 2:  
        #载体测序  
        #uha_dha测序
        sequencing_primer_df=uha_dha_sgRNA_df[[ 'Name', 'Region', 'UHA', 'UHA_size', 'DHA', 'DHA_size',
                                                'plasmid', 'promoter_N20_terminator',
                                                'promoter_N20_terminator_up', 'promoter_N20_terminator_down','seq_altered']]
        sequencing_primer_df['plasmid_sequencing_region'] = sequencing_primer_df['UHA'] + sequencing_primer_df['seq_altered'] + sequencing_primer_df['DHA'] + sequencing_primer_df['promoter_N20_terminator_up'] + sequencing_primer_df['promoter_N20_terminator'] 
        sequencing_primer_df['Region'] = sequencing_primer_df['Name']+';'+sequencing_primer_df['Region']
        sequencing_primer_template = sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
        plasmid_sequencing_primer_df,failture_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')

    temp = uha_dha_sgRNA_df[['Name','Region','UHA','DHA','promoter_N20_terminator','seq_altered']]
    temp['Region'] = temp['Name']+';'+temp['Region']
    temp = temp.drop(columns='Name')
    sequencing_primer_template = pd.merge(sequencing_primer_template,temp,on=['Region'],how='inner')

    return plasmid_sequencing_primer_df, sequencing_primer_template, failture_plasmid_sequencing_primer_df 

def two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df,promoter_seq,enzyme_df,enzyme_name):


    sgRNA_enzyme_df = enzyme_df[enzyme_df['name']==enzyme_name]
    sgRNA_enzyme_df = sgRNA_enzyme_df.reset_index(drop=True)
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
        sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df = two_plasmid_system_design_by_user_sgRNA_region(uha_dha_sgRNA_df,sgRNA_plasmid_seq,sgRNA_plasmid_region_seq,sgRNA_region_json,no_ccdb_plasmid,enzyme_df,enzyme_name)
    else:
        sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df = two_plasmid_system_design_by_sgRNA_no_user(no_ccdb_uha_dha_sgRNA_df, uha_dha_sgRNA_df, enzyme_df, enzyme_name)

    if ccdb_region_json != {}:        
        ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_ccdb_region(uha_dha_sgRNA_df,no_sgRNA_plasmid,ccdb_plasmid_seq,ccdb_plasmid_region_seq,ccdb_region_json, enzyme_df, enzyme_name)
    else:
        ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name)

    return  sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df, failture_sgRNA_plasmid_primer_df, failture_ccdb_plasmid_primer_df


def two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name):
        
        #执行用户什么都不指定，设计ccdb质粒引物   
        ccdb_primer_template_df = pd.DataFrame(columns=['plasmid_backbone_primer','plasmid_backbone'],data=[[f'ccdb_plasmid;primer',ccdB_plasmid_backbone]])
        ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = p_d_seq.design_primer(ccdb_primer_template_df,'plasmid_backbone_primer','plasmid_backbone','plasmid')
        if len(ccdb_plasmid_primer_df)>0:
            #ccdb质粒引物加接头
            ccdb_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(ccdb_plasmid_primer_df,enzyme_df,enzyme_name,stype='ccdb_plasmid_primer_joint')
            #提取引物的必要部分
            ccdb_plasmid_primer_df = ccdb_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]    

        return ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df

def two_plasmid_system_design_by_sgRNA_no_user(no_ccdb_uha_dha_sgRNA_df, uha_dha_sgRNA_df, enzyme_df, enzyme_name):
        
        #执行用户什么都不指定，设计sgRNA质粒引物
        sgRNA_primer_template_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','sgRNA_template','Rev Target sequence']]
        sgRNA_primer_template_df['Region'] = sgRNA_primer_template_df['Name'] +';'+ sgRNA_primer_template_df['Region']
        sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df  = p_d_seq.design_primer(sgRNA_primer_template_df,'Region','sgRNA_template','sgRNA')
        if len(sgRNA_plasmid_primer_df) > 0:

            uha_dha_sgRNA_df['Region'] = uha_dha_sgRNA_df['Name'] + ';' + uha_dha_sgRNA_df['Region']
            #sgRNA质粒引物加接头
            sgRNA_plasmid_primer_df = pd.merge(uha_dha_sgRNA_df[['Region','Target sequence','Rev Target sequence']],sgRNA_plasmid_primer_df,on=['Region'],how='inner')
            sgRNA_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(sgRNA_plasmid_primer_df,enzyme_df,enzyme_name,stype='sgRNA_plasmid_primer_joint')
            #提取关键信息
            sgRNA_plasmid_primer_df = sgRNA_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

        return sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df

def two_plasmid_system_design_by_user_sgRNA_region(uha_dha_sgRNA_df,sgRNA_plasmid_seq,sgRNA_plasmid_region_seq,sgRNA_region_json,no_ccdb_plasmid,enzyme_df,enzyme_name):
        
        #sgRNA载体质粒                   #启动子终止子若横跨零点有问题
        promoter_terminator_seq = sgRNA_plasmid_region_seq['promoter_seq'] + sgRNA_plasmid_region_seq['n20_coordinate_seq'] + sgRNA_plasmid_region_seq['terminator_seq']
        promoter_terminator_start = sgRNA_plasmid_seq.find(promoter_terminator_seq)
        promoter_terminator_end = promoter_terminator_start + len(promoter_terminator_seq)

        first_primer_position_in_promoter_terminator = promoter_terminator_seq.find(sgRNA_plasmid_region_seq['terminator_seq'])
        first_primer_start_position = promoter_terminator_start + first_primer_position_in_promoter_terminator
        sgRNA_plasmid_seq_len = len(sgRNA_plasmid_seq)

        sgRNA_distance_dict = region_2_distance(sgRNA_plasmid_seq_len, sgRNA_region_json,first_primer_start_position)
        state,sgRNA_plasmid_primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(first_primer_start_position, sgRNA_plasmid_seq, sgRNA_distance_dict, 20)
        
        sgRNA_plasmid_primer_df = pd.DataFrame()
        failture_sgRNA_plasmid_primer_df = pd.DataFrame()

        if state == 'success':
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
        elif state == 'failture':
            failture_sgRNA_plasmid_primer_df = pd.DataFrame(sgRNA_plasmid_primer_result_list)

        return sgRNA_plasmid_primer_df,failture_sgRNA_plasmid_primer_df

def two_plasmid_system_design_by_user_ccdb_region(uha_dha_sgRNA_df,no_sgRNA_plasmid,ccdb_plasmid_seq,ccdb_plasmid_region_seq,ccdb_region_json, enzyme_df, enzyme_name):

        #ccdb载体质粒
        ccdb_start = ccdb_plasmid_seq.find(ccdb_plasmid_region_seq['ccdb'])
        ccdb_end = ccdb_start + len(ccdb_plasmid_region_seq['ccdb'])
        first_primer_start_position = ccdb_end
        ccdb_plasmid_seq_len = len(ccdb_plasmid_seq)
        ccdb_distance_dict = region_2_distance(ccdb_plasmid_seq_len, ccdb_region_json,first_primer_start_position)

        state,ccdb_plasmid_primer_result_list = p_d_seq.design_primer_by_region_in_plasmid(first_primer_start_position,ccdb_plasmid_seq, ccdb_distance_dict,len(ccdb_plasmid_region_seq['ccdb']))
        
        failture_ccdb_plasmid_primer_df = pd.DataFrame()
        ccdb_plasmid_primer_df = pd.DataFrame()
        if state == 'success':
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
        elif state == 'failture':
            failture_ccdb_plasmid_primer_df = pd.DataFrame()

        return ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df


def two_plasmid_system_design_by_no_user(no_ccdb_uha_dha_sgRNA_df,ccdB_plasmid_backbone,enzyme_df,enzyme_name):
    #设计sgRNA质粒引物
    no_ccdb_primer_template_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','sgRNA_template','Rev Target sequence']]
    no_ccdb_primer_template_df['Region'] = no_ccdb_primer_template_df['Name'] + ';' + no_ccdb_primer_template_df['Region']
    sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df = p_d_seq.design_primer(no_ccdb_primer_template_df,'Region','sgRNA_template','sgRNA')
    no_ccdb_uha_dha_sgRNA_df['Region'] = no_ccdb_uha_dha_sgRNA_df['Name'] +';'+ no_ccdb_uha_dha_sgRNA_df['Region']
    
    #sgRNA质粒引物加接头
    if len(sgRNA_plasmid_primer_df)>0:
        sgRNA_plasmid_primer_df = pd.merge(no_ccdb_uha_dha_sgRNA_df[['Region','Target sequence','Rev Target sequence']],sgRNA_plasmid_primer_df,on=['Region'],how='inner')
        sgRNA_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(sgRNA_plasmid_primer_df,enzyme_df,enzyme_name,stype='sgRNA_plasmid_primer_joint')
    
    #设计ccdB质粒引物
    no_sgRNA_primer_template_df = pd.DataFrame(columns=['plasmid_backbone_primer','plasmid_backbone'],data=[[f'ccdb_plasmid;primer',ccdB_plasmid_backbone]])

    ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = p_d_seq.design_primer(no_sgRNA_primer_template_df,'plasmid_backbone_primer','plasmid_backbone','plasmid')
    
    #ccdb质粒引物加接头
    if len(ccdb_plasmid_primer_df)>0:
        ccdb_plasmid_primer_df = p_d_seq.add_joint_sgRNA_primer(ccdb_plasmid_primer_df,enzyme_df,enzyme_name,stype='ccdb_plasmid_primer_joint')
      
    return sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df, failture_sgRNA_plasmid_primer_df, failture_ccdb_plasmid_primer_df  

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
    state, sgRNA_primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         sgRNA_primer_json,
                                                         sgRNA_primer_position_json,
                                                         no_ccdb_plasmid,
                                                         sgRNA_plasmid_backbone,
                                                         n_20_label,
                                                         ccdb_label,
                                                         promoter_terminator_label)
                                                         
    sgRNA_plasmid_primer_df = pd.DataFrame()
    failture_sgRNA_plasmid_primer_df = pd.DataFrame()

    if state == 'success':
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
    elif state == 'failture':
        failture_sgRNA_plasmid_primer_df = sgRNA_primer_dict_df
    return sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df

def two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone, ccdb_primer_json, ccdb_primer_position_json, no_sgRNA_plasmid,enzyme_df,enzyme_name,n_20_label,ccdb_label,promoter_terminator_label):

    
    sgRNA_promoter_terminator_start = 0
    #对引物进行排序,确定所有引物
    state, ccdb_primer_dict_df, primers_sum = p_d_seq.sort_compose_primer(sgRNA_promoter_terminator_start,
                                                         ccdb_primer_json,
                                                         ccdb_primer_position_json,
                                                         no_sgRNA_plasmid,
                                                         ccdB_plasmid_backbone,
                                                         n_20_label,
                                                         ccdb_label,
                                                         promoter_terminator_label)
    ccdb_plasmid_primer_df = pd.DataFrame()
    failture_ccdb_plasmid_primer_df = pd.DataFrame()
    if state == 'success':
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
    elif state == 'failture':
        failture_ccdb_plasmid_primer_df = ccdb_primer_dict_df

    return ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df

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
        if sgRNA_failture_primer != {}:
            raise ValueError("primer error:", sgRNA_failture_primer)
    else:
        sgRNA_plasmid_primer_df,failture_sgRNA_plasmid_primer_df = two_plasmid_system_design_by_sgRNA_no_user(no_ccdb_uha_dha_sgRNA_df, uha_dha_sgRNA_df, enzyme_df, enzyme_name)
        
    if ccdb_primer_json !={}:  
        ccdb_primer_position_json, ccdb_failture_primer = p_d_seq.check_locate_primer(ccdB_plasmid_backbone, ccdb_primer_json)
        if ccdb_failture_primer != {}:
            raise ValueError("primer error:", ccdb_failture_primer)
    else:
        ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name)

    if  sgRNA_primer_position_json == {} and sgRNA_primer_json !={}:
        raise ValueError("primer error:", sgRNA_primer_json)
    elif sgRNA_primer_position_json !={} and sgRNA_primer_json != {}:
        #设计sgRNA primer
        sgRNA_plasmid_primer_df, failture_sgRNA_plasmid_primer_df = two_plasmid_system_design_by_user_primer_sgRNA( 
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
        raise ValueError("primer error:", ccdb_primer_json)
    elif ccdb_primer_position_json != {} and ccdb_primer_json != {}:
        #设计ccdb primer
        ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone,
                                                                                ccdb_primer_json,
                                                                                ccdb_primer_position_json,
                                                                                no_sgRNA_plasmid,
                                                                                enzyme_df,
                                                                                enzyme_name,
                                                                                n_20_label,
                                                                                ccdb_label,
                                                                                promoter_terminator_label)
                 
    return sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df, failture_sgRNA_plasmid_primer_df, failture_ccdb_plasmid_primer_df

def two_plasmid_system_sequencing_design_primer(no_ccdb_uha_dha_sgRNA_df, no_sgRNA_uha_dha_ccdb_df):

    #sgRNA质粒测序
    sgRNA_plasmid_sequencing_primer_df = no_ccdb_uha_dha_sgRNA_df[['Name','Region','plasmid','promoter_N20_terminator','promoter_N20_terminator_up','promoter_N20_terminator_down','seq_altered']]
    sgRNA_plasmid_sequencing_primer_df['plasmid_sequencing_region'] = sgRNA_plasmid_sequencing_primer_df['promoter_N20_terminator']
    sgRNA_plasmid_sequencing_primer_df['Region'] = sgRNA_plasmid_sequencing_primer_df['Name']+';'+sgRNA_plasmid_sequencing_primer_df['Region']
    sgRNA_plasmid_sequencing_primer_template = sgRNA_plasmid_sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
    sgRNA_plasmid_sequencing_primer_df,failture_sgRNA_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(sgRNA_plasmid_sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')
    #ccdb质粒载体测序
    ccdb_plasmid_sequencing_primer_df = no_sgRNA_uha_dha_ccdb_df[['Name','Region','UHA','DHA','seq_altered','plasmid']]
    ccdb_plasmid_sequencing_primer_df['plasmid_sequencing_region'] = ccdb_plasmid_sequencing_primer_df['UHA'] + ccdb_plasmid_sequencing_primer_df['seq_altered'] + ccdb_plasmid_sequencing_primer_df['DHA']
    ccdb_plasmid_sequencing_primer_df['Region'] =  ccdb_plasmid_sequencing_primer_df['Name'] + ';' + ccdb_plasmid_sequencing_primer_df['Region']
    ccdb_plasmid_sequencing_primer_template = ccdb_plasmid_sequencing_primer_df[['Name','Region','plasmid_sequencing_region','plasmid']]
    ccdb_plasmid_sequencing_primer_df, failture_ccdb_plasmid_sequencing_primer_df = p_d_seq.create_sequencing_primer(ccdb_plasmid_sequencing_primer_template,sr,'plasmid','plasmid_sequencing_region',seq_type='plasmid_seq')

    sgRNA_temp = no_ccdb_uha_dha_sgRNA_df[['Name','Region','promoter_N20_terminator']]
    sgRNA_temp['Region'] = sgRNA_temp['Name'] + ';' + sgRNA_temp['Region']
    ccdb_temp = no_sgRNA_uha_dha_ccdb_df[['Name','Region','UHA','DHA','seq_altered']]
    ccdb_temp['Region'] = ccdb_temp['Name'] + ';' + ccdb_temp['Region']

    sgRNA_plasmid_sequencing_primer_template = pd.merge(sgRNA_plasmid_sequencing_primer_template, sgRNA_temp, how='inner')
    ccdb_plasmid_sequencing_primer_template = pd.merge(ccdb_plasmid_sequencing_primer_template, ccdb_temp, how='inner')

    return sgRNA_plasmid_sequencing_primer_df, sgRNA_plasmid_sequencing_primer_template, ccdb_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_template, failture_sgRNA_plasmid_sequencing_primer_df,failture_ccdb_plasmid_sequencing_primer_df

#genome sequencing primer
def genome_sequencing_design_primer(info_input_df, uha_dha_df):

    #编辑基因组测序
    info_input_df1 = info_input_df[['name','region','seq_uha_max_whole','seq_dha_max_whole','uha_upstream','dha_downstream']].rename(columns={'name':'Name','region':'Region'})
    UHA_DHA_df = pd.merge(info_input_df1,uha_dha_df,how='inner')  #去除
    def work1(uha_upstream, seq_uha_max_whole, seq_altered, seq_dha_max_whole, dha_downstream):
        if seq_altered !='-':
            sequencing_template = uha_upstream + seq_uha_max_whole + seq_altered + seq_dha_max_whole + dha_downstream
        else:
            sequencing_template = uha_upstream + seq_uha_max_whole  + seq_dha_max_whole + dha_downstream
        return sequencing_template
    UHA_DHA_df['sequencing_template'] = UHA_DHA_df.apply(lambda x: work1(x['uha_upstream'], x['seq_uha_max_whole'], x['seq_altered'], x['seq_dha_max_whole'], x['dha_downstream']), axis=1)
    def work2(UHA, seq_altered, DHA):
        if seq_altered != '-':
            sequencing_region = UHA + seq_altered + DHA
        else:
            sequencing_region = UHA + DHA
        return sequencing_region
    UHA_DHA_df['sequencing_region'] = UHA_DHA_df.apply(lambda x: work2(x['UHA'], x['seq_altered'], x['DHA']), axis=1)
    ###
    UHA_DHA_df['Region'] =  UHA_DHA_df['Name']+';'+UHA_DHA_df['Region']
    genome_sequencing_primer_df,failture_genome_sequencing_primer_df = p_d_seq.create_sequencing_primer(UHA_DHA_df,sr,'sequencing_template','sequencing_region',seq_type='genome_seq')
    #测序模板
    genome_sequencing_template = UHA_DHA_df[['Region','sequencing_template','UHA','seq_altered','DHA']]

    return genome_sequencing_primer_df, genome_sequencing_template, failture_genome_sequencing_primer_df


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
                                output,
                                ccdb_label,
                                promoter_terminator_label,
                                n_20_label,
                                promoter_label,
                                failture_uha_dha_primer_df
                                ):
    
    gb = SeqIO.read(gb_path, "genbank")
    gb_name = gb.name
    
    #质粒引物的设计类型：1---用户指定范围，2----无需用户指定范围，3----用户指定额外引物  
    if plasmid_primer_desgin_type == 2:    
        #无需用户指定范围
        uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_p_df,failture_seq_altered_primer_df,type_kind,failture_n20up_primer_df, failture_n20down_primer_df = one_plasmid_system_pcr_design_primer(                                                                                                                            
                                                                                                                                        gb_path,
                                                                                                                                        info_df,
                                                                                                                                        uha_dha_sgRNA_df,
                                                                                                                                        uha_dha_info_primer_df,
                                                                                                                                        uha_dha_primer_df,
                                                                                                                                        enzyme_df,
                                                                                                                                        enzyme_name,
                                                                                                                                        plasmid_primer_desgin_type,
                                                                                                                                        sgRNA_region_seq_json,
                                                                                                                                        ccdb_label,
                                                                                                                                        promoter_terminator_label,
                                                                                                                                        n_20_label,
                                                                                                                                        promoter_label,
                                                                                                                                        )

        #修改 n20down_primer_p_df
        temp = n20down_primer_p_df.loc[:0]  
        temp.loc[0,'Region'] = 1
        n20down_primer_p_df = temp

    elif plasmid_primer_desgin_type == 1:
        #用户指定范围
        uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_p_df,failture_seq_altered_primer_df,type_kind,failture_n20up_primer_df, failture_n20down_primer_df = one_plasmid_system_pcr_design_primer(                                                                                                                            
                                                                                                                                        gb_path,
                                                                                                                                        info_df,
                                                                                                                                        uha_dha_sgRNA_df,
                                                                                                                                        uha_dha_info_primer_df,
                                                                                                                                        uha_dha_primer_df,
                                                                                                                                        enzyme_df,
                                                                                                                                        enzyme_name,
                                                                                                                                        plasmid_primer_desgin_type,
                                                                                                                                        sgRNA_region_seq_json,
                                                                                                                                        ccdb_label,
                                                                                                                                        promoter_terminator_label,
                                                                                                                                        promoter_label,
                                                                                                                                        n_20_label
                                                                                                                                        )
    elif plasmid_primer_desgin_type == 3:
         #用户指定额外引物  
        uha_dha_sgRNA_df,uha_primer_df,dha_primer_df,n20down_primer_p_df,n20up_primer_p_df, seq_altered_p_df,failture_seq_altered_primer_df,type_kind,failture_n20up_primer_df, failture_n20down_primer_df = one_plasmid_system_pcr_design_primer(                                                                                                                            
                                                                                                                                        gb_path,
                                                                                                                                        info_df,
                                                                                                                                        uha_dha_sgRNA_df,
                                                                                                                                        uha_dha_info_primer_df,
                                                                                                                                        uha_dha_primer_df,
                                                                                                                                        enzyme_df,
                                                                                                                                        enzyme_name,
                                                                                                                                        plasmid_primer_desgin_type,
                                                                                                                                        sgRNA_primer_json,
                                                                                                                                        ccdb_label,
                                                                                                                                        promoter_terminator_label,
                                                                                                                                        promoter_label,
                                                                                                                                        n_20_label
                                                                                                                                        )

    #设计质粒测序引物
    plasmid_sequencing_primer_df, sequencing_primer_template, failture_plasmid_sequencing_primer_df = one_plasmid_system_sequencing_design_primer(type_kind,uha_dha_sgRNA_df)
    
    #设计基因组测序引物
    genome_sequencing_primer_df, genome_sequencing_template, failture_genome_sequencing_primer_df = genome_sequencing_design_primer(info_input_df, uha_dha_df)

    #标准化，重命名
    if len(seq_altered_p_df)>0:
        df_common_list = su.rename_common_primer_df(n20up_primer_p_df, n20down_primer_p_df, seq_altered_p_df)
        n20up_primer_p_df, n20down_primer_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1], df_common_list[2]
        
    else:
        df_common_list = su.rename_common_primer_df(n20up_primer_p_df, n20down_primer_p_df)
        n20up_primer_p_df, n20down_primer_p_df = df_common_list[0], df_common_list[1]
        seq_altered_p_df = pd.DataFrame()
    
    uha_primer_df, dha_primer_df = su.rename_u_d_primer_df(uha_primer_df, dha_primer_df)
    df_sequencing_list = su.rename_sequencing_primer_df(plasmid_sequencing_primer_df, genome_sequencing_primer_df)
    plasmid_sequencing_primer_df, genome_sequencing_primer_df = df_sequencing_list[0], df_sequencing_list[1]


    #----------------------------生成gb文件用于可视化展示-------------------------------------------------------------------

   
    if len(n20up_primer_p_df) > 0 and len(n20down_primer_p_df) > 0:
        #筛选出：
        plasmid_primer_featrue_df = su.create_plasmid_primer_featrue_df(sequencing_primer_template,
                                                                    uha_primer_df,
                                                                    seq_altered_p_df,
                                                                    dha_primer_df,   
                                                                    n20up_primer_p_df)
        plasmid_primer_featrue_df = plasmid_primer_featrue_df.fillna('')
        joint_len, cut_seq_len = su.get_joint_by_enzyme(enzyme_df, enzyme_name)

        #为每个编辑区域创建gb文件
        pcr_gb_output = os.path.join(output,'one_plasmid_system_pcr_gb/')
        if not exists(pcr_gb_output):
            os.makedirs(pcr_gb_output)
        pcr_tsv_df = p_d_seq.create_gb_for_region(plasmid_primer_featrue_df, n20down_primer_p_df, joint_len, cut_seq_len, pcr_gb_output, type='sgRNA_ccdb')
    else:
        pcr_tsv_df = pd.DataFrame()


    #质粒测序引物模板
    plasmid_sequencing_template = sequencing_primer_template
    plasmid_sequencing_template.rename(columns={'Region':'ID','plasmid':'PLASMID','seq_altered':"SEQ_ALTERED",'promoter_N20_terminator':'PROMOTER_N20_TERMINATOR'},inplace=True)
    plasmid_seq_gb_output = os.path.join(output,'one_plasmid_system_plasmid_sequencing_gb/')
    if not exists(plasmid_seq_gb_output):
        os.makedirs(plasmid_seq_gb_output)
    
    #筛选出
    if len(plasmid_sequencing_template) > 0 and len(plasmid_sequencing_primer_df) > 0: 
        temp_df = pd.merge(plasmid_sequencing_template, plasmid_sequencing_primer_df,how='inner')
        plasmid_sequencing_template = temp_df[plasmid_sequencing_template.columns]
        plasmid_sequencing_primer_df = temp_df[plasmid_sequencing_primer_df.columns]
        plasmid_seq_tsv_df = p_d_seq.create_gb_for_sequencing_region(plasmid_sequencing_template, plasmid_sequencing_primer_df, plasmid_seq_gb_output, type='plasmid_sequencing')
    else:
        plasmid_seq_tsv_df = pd.DataFrame()

    #基因组测序引物模板
    genome_sequencing_template.rename(columns={'Region':'ID','sequencing_template':'PLASMID','seq_altered':"SEQ_ALTERED"},inplace=True)
    genome_seq_gb_output = os.path.join(output,'one_plasmid_system_genome_sequencing_gb/')
    if not exists(genome_seq_gb_output):
        os.makedirs(genome_seq_gb_output)

    genome_sequencing_primer_df_no_sequencing_target_length_df = genome_sequencing_primer_df[genome_sequencing_primer_df.columns[:-1]]
    
    #筛选
    if len(genome_sequencing_template) > 0 and len(genome_sequencing_primer_df_no_sequencing_target_length_df) > 0:
        temp_df = pd.merge(genome_sequencing_template, genome_sequencing_primer_df_no_sequencing_target_length_df,how='inner')
        genome_sequencing_template = temp_df[genome_sequencing_template.columns]
        genome_sequencing_primer_df_no_sequencing_target_length_df = temp_df[genome_sequencing_primer_df_no_sequencing_target_length_df.columns]
        genome_seq_tsv_df = p_d_seq.create_gb_for_sequencing_region(genome_sequencing_template, genome_sequencing_primer_df_no_sequencing_target_length_df, genome_seq_gb_output, type='genome_sequencing')
    else:
        genome_seq_tsv_df = pd.DataFrame()
    #合并三个df 
    if len(pcr_tsv_df)>0 and len(plasmid_seq_tsv_df)>0 and len(genome_seq_tsv_df)>0:
        tsv_df = pd.merge(pcr_tsv_df,plasmid_seq_tsv_df,on='name',how='inner').merge(genome_seq_tsv_df,on='name',how='inner')
    else:
        tsv_df = pd.DataFrame()
    tsv_df.to_csv(os.path.join(output,'one_plasmid_system_gb_visualization.tsv'), index=False, sep='\t') 
    #-----------------------------生成gb文件用于可视化展示------------------------------------------------------------------------------------------   

    parent_output  = output
    output = output+'/one_plasmid_system_result/'
    if not exists(output):
        os.makedirs(output)

    xlsx_file = os.path.join(
        output,
        'one_plasmid_design_result.xlsx' 
    )
    order_file = os.path.join(
        output,
        'one_plasmid_primer_order.xlsx'
    )
    xlsx_F_file = os.path.join(
        output,
        'one_plasmid_design_F_Result.xlsx' 
    )

    #0：输出引物合成订单
    unique_all_primer_df = order.merge_primer(uha_primer_df, dha_primer_df, n20up_primer_p_df, n20down_primer_p_df, seq_altered_p_df, plasmid_sequencing_primer_df, genome_sequencing_primer_df)
    order.create_orders(unique_all_primer_df, orders_path = order_file)


    #1: 生成质粒合成订单
    plasmid_order_df = create_plasmid_order(uha_dha_sgRNA_df, name='', plasmid_type='one_plasmid_system')


    #补充:
    if len(uha_dha_sgRNA_df) > 0 and len(genome_sequencing_primer_df) > 0:
        genome_sequencing_primer_df = add_WT_genome_sequencing_primer_df(uha_dha_sgRNA_df, genome_sequencing_primer_df)

    #2：输出成功引物excel文件
    with pd.ExcelWriter(xlsx_file) as writer:  
        uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
        dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
        seq_altered_p_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')
        n20up_primer_p_df.to_excel(writer,sheet_name = 'Primer_sgRNA_fragment',index_label='No.')  
        n20down_primer_p_df.to_excel(writer,sheet_name = 'Primer_plasmid_backbone',index_label='No.')
        plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
        genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')
        # plasmid_order_df.to_excel(writer,sheet_name = 'Plasmid_synthesis_order',index_label='No.')
    
    #3.输出失败引物excel文件
    #分离uha，dha
    if len(failture_uha_dha_primer_df) >0 :
        failture_uha_primer_df = failture_uha_dha_primer_df[failture_uha_dha_primer_df['type']=='uha']
        failture_uha_primer_df = failture_uha_primer_df.drop(columns='type')
        failture_uha_primer_df.reset_index(drop=True, inplace=True) 
        failture_dha_primer_df = failture_uha_dha_primer_df[failture_uha_dha_primer_df['type']=='dha']
        failture_dha_primer_df = failture_dha_primer_df.drop(columns='type')
        failture_dha_primer_df.reset_index(drop=True, inplace=True)
    else:   
        failture_uha_primer_df = pd.DataFrame()
        failture_dha_primer_df = pd.DataFrame()

    failture_plasmid_sequencing_primer_df.reset_index(drop=True, inplace=True)
    failture_genome_sequencing_primer_df.reset_index(drop=True, inplace=True)
    


    with pd.ExcelWriter(xlsx_F_file) as writer: 

        failture_uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
        failture_dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
        failture_seq_altered_primer_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')
        failture_n20up_primer_df.to_excel(writer,sheet_name = 'Primer_sgRNA_fragment',index_label='No.')
        failture_n20down_primer_df.to_excel(writer,sheet_name = 'Primer_plasmid_backbone',index_label='No.')
        failture_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
        failture_genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')   


    #复制one_plasmid_system_pcr_gb到one_plasmid_system_result中
    os.system(f"cp -r {os.path.join(parent_output,'one_plasmid_system_pcr_gb')} {output}")
    os.system(f"cp -r {os.path.join(parent_output,'blast/Evaluation_result.txt')} {output}") 

    #输出一个zip文件夹     
    su.zip_ya(output, os.path.join(parent_output,'one_plasmid_system_result.zip'),num=0)

    return os.path.join(parent_output,'one_plasmid_system_result.zip')                                              

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
                                 promoter_label,
                                 output,
                                 failture_uha_dha_primer_df
                                 ):                                     

    no_ccdb_uha_dha_sgRNA_df,promoter_seq, sgRNA_plasmid_backbone, promoter_seq, terminator_seq, sgRNA_promoter_terminator = p_d_seq.create_new_plasmid(no_ccdb_plasmid, uha_dha_sgRNA_df.copy(), ccdb_label, promoter_terminator_label, n_20_label, promoter_label)
    no_sgRNA_uha_dha_ccdb_df, ccdB_plasmid_backbone, ccdB_promoter_terminator_up_seq = p_d_seq.create_new_plasmid(no_sgRNA_plasmid, uha_dha_sgRNA_df.copy(), ccdb_label, promoter_terminator_label, n_20_label)
    no_ccdb_uha_dha_sgRNA_df['promoter_seq'] = promoter_seq
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
        sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df, failture_sgRNA_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_region(
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
            sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df, failture_sgRNA_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_no_user(no_ccdb_uha_dha_sgRNA_df, ccdB_plasmid_backbone,enzyme_df,enzyme_name)
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
        sgRNA_plasmid_primer_df, ccdb_plasmid_primer_df, failture_sgRNA_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer(
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
            ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_ccdb_region(uha_dha_sgRNA_df,no_sgRNA_plasmid,ccdb_plasmid_seq,ccdb_plasmid_region_seq,ccdb_region_json, enzyme_df, enzyme_name)

    elif    plasmid_primer_desgin_type == 2 and method == 'OLIGO':
            #设计sgRNA质粒
            enzymeCutSeq_and_N20_df = two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df, promoter_seq, enzyme_df, enzyme_name)
            #设计ccdb质粒
            ccdb_plasmid_primer_df, failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_ccdb_no_user(ccdB_plasmid_backbone,enzyme_df,enzyme_name)

    elif    plasmid_primer_desgin_type == 3 and method == 'OLIGO':
            #设计sgRNA质粒
            enzymeCutSeq_and_N20_df = two_plasmid_system_n20_enzyme_cut_seq(no_ccdb_uha_dha_sgRNA_df, promoter_seq, enzyme_df, enzyme_name)
            #设计ccdb质粒
            ccdb_primer_position_json, ccdb_failture_primer = p_d_seq.check_locate_primer(ccdB_plasmid_backbone, ccdb_primer_json)

            if ccdb_primer_position_json == {} and ccdb_primer_json != {}:
                    return ccdb_failture_primer
            
            elif    ccdb_primer_position_json != {} and ccdb_primer_json != {}:
                    #设计ccdb primer
                    ccdb_plasmid_primer_df,failture_ccdb_plasmid_primer_df = two_plasmid_system_design_by_user_primer_ccdb(ccdB_plasmid_backbone,
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
    if len(seq_altered_primer_template) > 0:
        seq_altered_primer_template['Region'] = seq_altered_primer_template['Name'] +';'+ seq_altered_primer_template['Region']
        seq_altered_primer_df, failture_seq_altered_primer_df = p_d_seq.design_primer(seq_altered_primer_template,'Region','seq_altered',stype='seq_altered')
        #seq_altered_primer加接头
        seq_altered_primer_df = p_d_seq.add_joint_sgRNA_primer(seq_altered_primer_df, enzyme_df, enzyme_name, stype='seq_altered_primer_joint')

        #提取变化序列引物
        seq_altered_p_df = seq_altered_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

    else:
        seq_altered_primer_df = pd.DataFrame()
        seq_altered_p_df = seq_altered_primer_df
        failture_seq_altered_primer_df = pd.DataFrame()

    #给引物加接头
    #给同源臂引物加接头
    if len(uha_dha_info_primer_df) > 0:
        uha_dha_primer_df = p_d_seq.add_joint_sgRNA_primer(uha_dha_info_primer_df, enzyme_df, enzyme_name, ccdB_plasmid_backbone, stype = 'u_d_primer_joint')
        #分别提取，修饰后的uha、dha引物
        in_col = ['Name','Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint",'product_value_joint','product_size_joint']
        ou_col = ['Name','Region',"u_primer_f_seq_(5'-3')","u_primer_r_seq_(5'-3')",'UHA','UHA_size',"d_primer_f_seq_(5'-3')","d_primer_r_seq_(5'-3')",'DHA','DHA_size']
        uha_primer_df, dha_primer_df = p_d_seq.create_uha_dha_primer_df(uha_dha_primer_df, in_col, ou_col)
        uha_primer_df, dha_primer_df = su.rename_u_d_primer_df(uha_primer_df,dha_primer_df)

    #提取引物的必要部分
    if method == 'PCR' and len(sgRNA_plasmid_primer_df) > 0:
        sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]

    elif  method == 'PCR' and len(sgRNA_plasmid_primer_df) == 0:
        sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df

    
   
    if len(ccdb_plasmid_primer_df) > 0:
        ccdb_plasmid_p_df = ccdb_plasmid_primer_df[['Region',r"primer_f_seq_(5'-3')_joint",r"primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
    elif len(ccdb_plasmid_primer_df) == 0:
        ccdb_plasmid_p_df = ccdb_plasmid_primer_df


    #设计质粒测序引物  
    sgRNA_plasmid_sequencing_primer_df, sgRNA_plasmid_sequencing_primer_template, ccdb_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_template, failture_sgRNA_plasmid_sequencing_primer_df,failture_ccdb_plasmid_sequencing_primer_df = two_plasmid_system_sequencing_design_primer(no_ccdb_uha_dha_sgRNA_df,no_sgRNA_uha_dha_ccdb_df)

    #设计基因组测序引物
    genome_sequencing_primer_df, genome_sequencing_template, failture_genome_sequencing_primer_df = genome_sequencing_design_primer(info_input_df, uha_dha_df)

    #标准化，重命名
    if method == 'PCR':
        if len(sgRNA_plasmid_p_df) > 0 and len(ccdb_plasmid_p_df)>0 and len(seq_altered_p_df)>0:
            df_common_list = su.rename_common_primer_df(sgRNA_plasmid_p_df,ccdb_plasmid_p_df,seq_altered_p_df)
            sgRNA_plasmid_p_df, ccdb_plasmid_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1], df_common_list[2]

        elif len(sgRNA_plasmid_p_df) > 0 and len(ccdb_plasmid_p_df)>0 and len(seq_altered_p_df)==0:
            df_common_list = su.rename_common_primer_df(sgRNA_plasmid_p_df, ccdb_plasmid_p_df)
            sgRNA_plasmid_p_df, ccdb_plasmid_p_df = df_common_list[0], df_common_list[1]  
            seq_altered_p_df = seq_altered_primer_df

        elif len(sgRNA_plasmid_p_df) > 0 and len(ccdb_plasmid_p_df)==0 and len(seq_altered_p_df)==0:
            df_common_list = su.rename_common_primer_df(sgRNA_plasmid_p_df)
            sgRNA_plasmid_p_df = df_common_list[0]
            seq_altered_p_df = seq_altered_primer_df
            ccdb_plasmid_p_df = ccdb_plasmid_primer_df
        elif len(ccdb_plasmid_p_df) > 0 and len(sgRNA_plasmid_p_df)==0 and len(seq_altered_p_df)==0:
            df_common_list = su.rename_common_primer_df(ccdb_plasmid_p_df)
            ccdb_plasmid_p_df = df_common_list[0]
            seq_altered_p_df = seq_altered_primer_df
            sgRNA_plasmid_p_df = sgRNA_plasmid_primer_df
    else:
        if len(ccdb_plasmid_p_df) >0 and len(seq_altered_p_df)>0:
            df_common_list = su.rename_common_primer_df(ccdb_plasmid_p_df,seq_altered_p_df)
            ccdb_plasmid_p_df, seq_altered_p_df = df_common_list[0], df_common_list[1]
        elif len(ccdb_plasmid_p_df)>0 and len(seq_altered_p_df)==0:
            df_common_list = su.rename_common_primer_df(ccdb_plasmid_p_df)
            ccdb_plasmid_p_df = df_common_list[0]
            seq_altered_p_df = seq_altered_primer_df
        elif len(ccdb_plasmid_p_df) == 0 and len(seq_altered_p_df) == 0:
            ccdb_plasmid_p_df = ccdb_plasmid_primer_df
            seq_altered_p_df = seq_altered_primer_df
        elif len(ccdb_plasmid_p_df) == 0 and len(seq_altered_p_df) > 0:
            df_common_list = su.rename_common_primer_df(seq_altered_p_df)
            seq_altered_p_df = df_common_list[0]
            ccdb_plasmid_p_df = ccdb_plasmid_primer_df


    if len(sgRNA_plasmid_sequencing_primer_df)>0 and len(ccdb_plasmid_sequencing_primer_df)>0 and len(genome_sequencing_primer_df)>0:
        df_sequencing_list = su.rename_sequencing_primer_df(sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df, genome_sequencing_primer_df)
        sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df, genome_sequencing_primer_df = df_sequencing_list[0], df_sequencing_list[1], df_sequencing_list[2]
    elif len(sgRNA_plasmid_sequencing_primer_df)>0 and len(ccdb_plasmid_sequencing_primer_df) >0 and len(genome_sequencing_primer_df)==0:
        df_sequencing_list = su.rename_sequencing_primer_df(sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df)
        sgRNA_plasmid_sequencing_primer_df, ccdb_plasmid_sequencing_primer_df = df_sequencing_list[0], df_sequencing_list[1]
    elif len(sgRNA_plasmid_sequencing_primer_df)>0 and len(ccdb_plasmid_sequencing_primer_df) ==0 and len(genome_sequencing_primer_df)==0:
        df_sequencing_list = su.rename_sequencing_primer_df(sgRNA_plasmid_sequencing_primer_df)
        sgRNA_plasmid_sequencing_primer_df = df_sequencing_list[0]
    elif len(sgRNA_plasmid_sequencing_primer_df)==0 and len(ccdb_plasmid_sequencing_primer_df) >0 and len(genome_sequencing_primer_df)==0:
        df_sequencing_list = su.rename_sequencing_primer_df(ccdb_plasmid_sequencing_primer_df)
        ccdb_plasmid_sequencing_primer_df = df_sequencing_list[0]
    elif len(sgRNA_plasmid_sequencing_primer_df)==0 and len(ccdb_plasmid_sequencing_primer_df) ==0 and len(genome_sequencing_primer_df)>0:
        df_sequencing_list = su.rename_sequencing_primer_df(genome_sequencing_primer_df)
        genome_sequencing_primer_df = df_sequencing_list[0]
    
    
    #------------------------------------------生成gb文件用于引物的可视化展示-------------------------------------------------------------
    #--------------------------PCR----------------------生成sgRNA_gb文件-------------------------------------------------------------------

    if len(sgRNA_plasmid_sequencing_primer_template) > 0:
         #为每个编辑区域创建gb文件
        gb_output = os.path.join(output,'two_plasmid_system_pcr_gb/')
        if not exists(gb_output):
            os.makedirs(gb_output)   

        joint_len, cut_seq_len = su.get_joint_by_enzyme(enzyme_df,enzyme_name)

        sgRNA_plasmid_sequencing_primer_template = sgRNA_plasmid_sequencing_primer_template[['Region','plasmid','promoter_N20_terminator']].rename(columns={'Region':'ID','plasmid':"PLASMID","promoter_N20_terminator":"PROMOTER_N20_TERMINATOR"})
        if method == 'PCR':
            if len(sgRNA_plasmid_p_df) >0 :
                sgRNA_plasmid_primer = sgRNA_plasmid_p_df[['ID', 'PRIMER_LEFT_WHOLE_SEQUENCE', 'PRIMER_RIGHT_WHOLE_SEQUENCE']]
                type = 'sgRNA'
                sgRNA_pcr_tsv_df = p_d_seq.create_gb_for_region(sgRNA_plasmid_sequencing_primer_template, sgRNA_plasmid_primer, joint_len, cut_seq_len, gb_output,type)
        elif method == 'OLIGO':
            sgRNA_plasmid_primer = enzymeCutSeq_and_N20_df[['ID','Target sequence']]
            type = 'enzyme_cut'
            sgRNA_pcr_tsv_df = p_d_seq.create_gb_for_region(sgRNA_plasmid_sequencing_primer_template, sgRNA_plasmid_primer, joint_len, cut_seq_len, gb_output,type)
        
       
        # if method == 'PCR':
        #     type = 'sgRNA'
        #     if len(sgRNA_plasmid_primer)>0:
        #         sgRNA_pcr_tsv_df = p_d_seq.create_gb_for_region(sgRNA_plasmid_sequencing_primer_template, sgRNA_plasmid_primer, joint_len, cut_seq_len, gb_output,type)
        # elif method == 'OLIGO':
        #     type = 'enzyme_cut'
        #     sgRNA_pcr_tsv_df = p_d_seq.create_gb_for_region(sgRNA_plasmid_sequencing_primer_template, sgRNA_plasmid_primer, joint_len, cut_seq_len, gb_output,type)
    else:
        sgRNA_pcr_tsv_df = pd.DataFrame(columns=['name'])


    #---------------------------PCR---------------------生成ccdb_gb文件---------------------------------------------------------------------
    # plasmid_primer_featrue_df = ccdb_plasmid_sequencing_primer_template[['Region','plasmid']].rename(columns={'Region':'ID','plasmid':"PLASMID"})
    # [['Name','Region','UHA','DHA','seq_altered']]
    if len(uha_primer_df) >0 and len(dha_primer_df)>0 and len(ccdb_plasmid_p_df)>0:   
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
        ccdb_pcr_tsv_df = p_d_seq.create_gb_for_region(plasmid_primer_featrue_df, ccdb_plasmid_p_df, joint_len, cut_seq_len, gb_output,type='ccdb')
    else:
        ccdb_pcr_tsv_df = pd.DataFrame(columns=['name'])


    #---------------------plasmid-------SEQUENCING------------------------------------------------------------------------------------------------------------
    #为每个编辑区域创建gb文件
    gb_output = os.path.join(output,'two_plasmid_system_plasmid_sequencing_gb/')
    if not exists(gb_output):
        os.makedirs(gb_output)

    sgRNA_plasmid_primer_featrue_df = sgRNA_plasmid_sequencing_primer_template
    if len(sgRNA_plasmid_primer_featrue_df) > 0 and len(sgRNA_plasmid_sequencing_primer_df) > 0:
        sgRNA_plasmid_primer_featrue_df.rename(columns={'Region':"ID", 'plasmid':"PLASMID", "promoter_N20_terminator":"PROMOTER_N20_TERMINATOR"},inplace=True)
        temp_df = pd.merge(sgRNA_plasmid_primer_featrue_df, sgRNA_plasmid_sequencing_primer_df,how='inner')
        sgRNA_plasmid_primer_featrue_df = temp_df[sgRNA_plasmid_primer_featrue_df.columns]  
        sgRNA_plasmid_sequencing_primer_df = temp_df[sgRNA_plasmid_sequencing_primer_df.columns]
        sgRNA_plasmid_sequencing_tsv_df = p_d_seq.create_gb_for_sequencing_region(sgRNA_plasmid_primer_featrue_df, sgRNA_plasmid_sequencing_primer_df, gb_output, type='sgRNA_plasmid_sequencing')
    else:
        sgRNA_plasmid_sequencing_tsv_df = pd.DataFrame(columns=['name'])


    if len(ccdb_plasmid_sequencing_primer_template) and len(ccdb_plasmid_sequencing_primer_df) > 0:
        ccdb_plasmid_primer_featrue_df = ccdb_plasmid_sequencing_primer_template[['Region','plasmid','UHA','DHA']].rename(columns={'Region':'ID','plasmid':"PLASMID"})
        temp_df = pd.merge(ccdb_plasmid_primer_featrue_df, ccdb_plasmid_sequencing_primer_df,how='inner')
        ccdb_plasmid_primer_featrue_df = temp_df[ccdb_plasmid_primer_featrue_df.columns]
        ccdb_plasmid_sequencing_primer_df = temp_df[ccdb_plasmid_sequencing_primer_df.columns]
        ccdb_plasmid_sequencing_tsv_df = p_d_seq.create_gb_for_sequencing_region(ccdb_plasmid_primer_featrue_df, ccdb_plasmid_sequencing_primer_df, gb_output, type='ccdb_plasmid_sequencing')
    else:
        ccdb_plasmid_sequencing_tsv_df = pd.DataFrame(columns=['name'])
    
   
    #-------------------genome----------SEQUENCING------------------------------------------------------------------
    gb_output = os.path.join(output,'two_plasmid_system_genome_sequencing_gb/')
    if not exists(gb_output):
        os.makedirs(gb_output)

    if len(genome_sequencing_template) >0 and len(genome_sequencing_primer_df) >0:
        genome_sequencing_template.rename(columns={'Region':'ID','sequencing_template':'PLASMID','seq_altered':"SEQ_ALTERED"},inplace=True)
        temp_df = pd.merge(genome_sequencing_template, genome_sequencing_primer_df,how='inner')
        genome_sequencing_template = temp_df[ genome_sequencing_template.columns ]
        genome_sequencing_primer_df = temp_df[ genome_sequencing_primer_df.columns ]
        genome_sequencing_tsv_df = p_d_seq.create_gb_for_sequencing_region(genome_sequencing_template, genome_sequencing_primer_df, gb_output, type='genome_sequencing')
    else:
        genome_sequencing_tsv_df = pd.DataFrame(columns=['name'])


    #生成tsv
    # print(sgRNA_tsv_df,'\n',ccdb_tsv_df)  
    pcr_df = pd.merge(sgRNA_pcr_tsv_df, ccdb_pcr_tsv_df, on='name')  
    sequencing_df = pd.merge(sgRNA_plasmid_sequencing_tsv_df, ccdb_plasmid_sequencing_tsv_df,on='name')
    tsv_df = pcr_df.merge(sequencing_df,on='name').merge(genome_sequencing_tsv_df,on='name')
    tsv_df.to_csv(os.path.join(output,'two_plasmid_system_gb_visualization.tsv'), index=False, sep='\t')
    #-------------------------------------------------------------------------------------------------------------------------------------------


    #1: 生成质粒合成订单
    ccdb_plasmid_order_df = create_plasmid_order(no_sgRNA_uha_dha_ccdb_df, name='ccdb_plasmid', plasmid_type='ccdb')
    sgRNA_plasmid_order_df = create_plasmid_order(no_ccdb_uha_dha_sgRNA_df, name='sgRNA_plasmid', plasmid_type='sgRNA')

    parent_output  = output   
    output = output+'/two_plasmid_system_result/'
    if not exists(output):
        os.makedirs(output)

    #输出引物订单
    order_file = os.path.join(
        output,
        'two_plasmid_primer_order.xlsx'
    )
    if method == "PCR":
        unique_all_primer_df = order.merge_primer(uha_primer_df,dha_primer_df,sgRNA_plasmid_p_df,ccdb_plasmid_p_df,seq_altered_p_df,sgRNA_plasmid_sequencing_primer_df,ccdb_plasmid_sequencing_primer_df,genome_sequencing_primer_df)
        order.create_orders(unique_all_primer_df, orders_path = order_file)
    else:
        unique_all_primer_df = order.merge_primer(uha_primer_df,dha_primer_df,enzymeCutSeq_and_N20_df,ccdb_plasmid_p_df,seq_altered_p_df,sgRNA_plasmid_sequencing_primer_df,ccdb_plasmid_sequencing_primer_df,genome_sequencing_primer_df)
        order.create_orders(unique_all_primer_df, orders_path = order_file)

    #输出引物     
    xlsx_file = os.path.join(
        output,
        'two_plasmid_design_result.xlsx'
    )
    xlsx_F_file = os.path.join(
        output,
        'two_plasmid_design_F_Result.xlsx' 
    )


     #3.输出失败引物excel文件
    #分离uha，dha
    if len(failture_uha_dha_primer_df) > 0:
        failture_uha_primer_df = failture_uha_dha_primer_df[failture_uha_dha_primer_df['type']=='uha']
        failture_uha_primer_df = failture_uha_primer_df.drop(columns='type')
        failture_uha_primer_df.reset_index(drop=True, inplace=True)

        failture_dha_primer_df = failture_uha_dha_primer_df[failture_uha_dha_primer_df['type']=='dha']
        failture_dha_primer_df = failture_dha_primer_df.drop(columns='type')
        failture_dha_primer_df.reset_index(drop=True, inplace=True)
    else:
        failture_uha_primer_df = pd.DataFrame()
        failture_dha_primer_df = pd.DataFrame()

    failture_sgRNA_plasmid_sequencing_primer_df.reset_index(drop=True, inplace=True)
    failture_ccdb_plasmid_sequencing_primer_df.reset_index(drop=True, inplace=True)
    
    if method == 'PCR':
        #成功
        with pd.ExcelWriter(xlsx_file) as writer:  
            uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            seq_altered_p_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')
            sgRNA_plasmid_p_df.to_excel(writer,sheet_name = 'Primer_sgRNA_fragment',index_label='No.')  
            ccdb_plasmid_p_df.to_excel(writer,sheet_name = 'Primer_plasmid_backbone',index_label='No.')
           
            sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
            ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
            genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')

            # ccdb_plasmid_order_df.to_excel(writer,sheet_name = 'HR_Plasmid_synthesis_order',index_label='No.')
            # sgRNA_plasmid_order_df.to_excel(writer,sheet_name = 'SgRNA_Plasmid_synthesis_order',index_label='No.') 

        #失败
        with pd.ExcelWriter(xlsx_F_file) as writer: 
            failture_uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            failture_dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            failture_seq_altered_primer_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')

            failture_sgRNA_plasmid_primer_df.to_excel(writer,sheet_name = 'Primer_sgRNA_fragment',index_label='No.')
            failture_ccdb_plasmid_primer_df.to_excel(writer,sheet_name = 'Primer_plasmid_backbone',index_label='No.')

            failture_sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
            failture_ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
            failture_genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.') 
    else:
        #成功
        with pd.ExcelWriter(xlsx_file) as writer:  
            uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            seq_altered_p_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')
            enzymeCutSeq_and_N20_df.to_excel(writer,sheet_name = 'Primer_sgRNA_fragment',index_label='No.')  
            ccdb_plasmid_p_df.to_excel(writer,sheet_name = 'Primer_plasmid_backbone',index_label='No.')  
            sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
            ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
            genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')

        
         #失败
        with pd.ExcelWriter(xlsx_F_file) as writer: 
            failture_uha_primer_df.to_excel(writer,sheet_name = 'Primer_UHA',index_label='No.')
            failture_dha_primer_df.to_excel(writer,sheet_name = 'Primer_DHA',index_label='No.')
            failture_seq_altered_primer_df.to_excel(writer,sheet_name = 'Primer_inserted_fragment',index_label='No.')

            failture_sgRNA_plasmid_primer_df.to_excel(writer,sheet_name = 'Primer_sgRNA_fragment',index_label='No.')
            failture_ccdb_plasmid_primer_df.to_excel(writer,sheet_name = 'Primer_plasmid_backbone',index_label='No.')

            failture_sgRNA_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P1',index_label='No.')
            failture_ccdb_plasmid_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_P2',index_label='No.')
            failture_genome_sequencing_primer_df.to_excel(writer,sheet_name = 'Test_primer_G',index_label='No.')  

        




    #复制one_plasmid_system_pcr_gb到one_plasmid_system_result中
    os.system(f"cp -r {os.path.join(parent_output,'two_plasmid_system_pcr_gb')} {output}")
    os.system(f"cp -r {os.path.join(parent_output,'blast/Evaluation_result.txt')} {output}")

    su.zip_ya(output, os.path.join(parent_output,'two_plasmid_system_result.zip'),num=0)

    return os.path.join(parent_output,'two_plasmid_system_result.zip') 


def read_chopchopInput_add_uha_dha(genome_path,chopchop_input,uha_dha_params):
    max_left_arm_seq_length = uha_dha_params['max_left_arm_seq_length']
    max_right_arm_seq_length = uha_dha_params['max_right_arm_seq_length']

    info_input_df = su.del_Unnamed(pd.read_csv(chopchop_input))
    info_input_df.columns = [i.lower() for i in info_input_df.columns]


    def work(mun_id, geneid, mutation_pos_index, ref):


        gene_start = mutation_pos_index
        gene_end = mutation_pos_index + len(ref)
        print(gene_start,gene_end)

        
        if mutation_pos_index - max_left_arm_seq_length < 0:
            error_message = "The length of upstream sequence of manipulation site of " + mun_id + " must be larger than sum of 'Max Length of UHA' and 'Max Length of UIS'."
            return error_message,error_message,error_message,error_message

        record = su.extract_seq_from_genome(genome_path,geneid)

        seq_uha_max_whole = str(record[
                        gene_start - max_left_arm_seq_length : gene_start 
                        ])
        seq_dha_max_whole = str(record[
                                gene_end : gene_end + max_right_arm_seq_length
                                ])


        uha_upstream = str(  
                        record[
                            gene_start - max_left_arm_seq_length - 200 : gene_start - max_left_arm_seq_length
                        ]
                    )
        dha_downstream=str(
                        record[
                            gene_end + max_right_arm_seq_length  : gene_end + max_right_arm_seq_length  + 200
                        ]
                    )
        return  uha_upstream, dha_downstream, seq_uha_max_whole, seq_dha_max_whole

    info_df = su.lambda2cols(info_input_df,lambdaf=work,in_coln=['name','geneid','mutation_pos_index','ref'],to_colns=['uha_upstream','dha_downstream','seq_uha_max_whole','seq_dha_max_whole'])

    return info_df    


def check_quality_control(plasmid_backbone,seq_json):
    
    failture_seq_json = {}
   
    for k,v in seq_json.items():
        plasmid_backbone = plasmid_backbone.upper()


        start = plasmid_backbone.find(v.upper())
        if start == -1:
            start = plasmid_backbone.find(su.revComp(v.upper()))

        errorMessage = ""


         #是否是DNA序列
        if not su.is_dna(v):
            errorMessage = f"{v}: The sequence is not a DNA sequence" 
            failture_seq_json.update({k:errorMessage})
            raise ValueError(errorMessage)  

        if start == -1:
            failture_seq_json.update({k: f'{v}:The sequence is not on the plasmid, please enter the sequence on the sense chain of the plasmid'})
            errorMessage = f'{v}: The sequence is not on the plasmid, please enter the sequence on the sense chain of the plasmid'
            failture_seq_json.update({k:errorMessage})
            raise ValueError(errorMessage)  
        else:
            #判断 5 < 序列长度 < 300
            if len(v) < 5  or  len(v) > 300:
                errorMessage = f"{v}: The sequence length should be between 5 and 300" 
                failture_seq_json.update({k:errorMessage}) 
                raise ValueError(errorMessage)  
           
            # 去取质粒骨架，引物的位置距离质粒骨架起始的距离大于50bp，引物的位置距离质粒骨架终止的距离大于50bp
            if start < 50 or len(plasmid_backbone) - (start + len(v)) < 50:
                errorMessage = f'{v} The distance between the position of the sequence and the start of the plasmid skeleton should be greater than 50 bp, and the distance between the position of the sequence and the end of the plasmid skeleton should be greater than 50 bp'
                failture_seq_json.update({k:errorMessage})
                raise ValueError(errorMessage)  
        

    return failture_seq_json

  
def check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label='gRNA'):
    gb = SeqIO.read(gb_path, "genbank") 

    selected_feature_labels = []
    # selected_feature_type = 'misc_feature'
    # selected_feature_key = 'gRNA'

    for feature in gb.features:
        if feature.type == selected_feature_type and 'label' in feature.qualifiers:
            label_value = feature.qualifiers['label'][0]
            if target_gene_label == label_value:
                selected_feature_labels.append(label_value)

    if len(selected_feature_labels) > 1:
        raise ValueError(f"{target_gene_label}:There is duplication in your plasmid feature labels")
    elif len(selected_feature_labels) == 0:
        raise ValueError(f"{target_gene_label}:Your plasmid feature label does not exist")  


def check_plasmid(gb_path, ccdb_label='', promoter_terminator_label='', n_20_label=''):

    try:
        gb = SeqIO.read(gb_path, "genbank")
    except Exception as e:
        raise ValueError(f"There is a problem with your plasmid gb file and an error occurred during the parsing process:{e}")
    gb_seq = str(gb.seq)
    #get coordinate
    ccdb_coordinate = su.get_feature_coordinate(ccdb_label,gb_path,selected_feature_type='misc_feature') 
    #N20
    n20_coordinate = su.get_feature_coordinate(n_20_label,gb_path,selected_feature_type='misc_feature')

    # #gRNA
    # gRNA_coordinate = su.get_feature_coordinate(promoter_terminator_label,gb_path)

    if ccdb_coordinate != (-1,-1) and n20_coordinate !=(-1,-1):
        check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label=ccdb_label)
        check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label=promoter_terminator_label)
        check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label=n_20_label)
        return 'one_plasmid_file_path'
    elif ccdb_coordinate==(-1,-1) and n20_coordinate !=(-1,-1):
        check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label=promoter_terminator_label)
        check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label=n_20_label)
        return 'no_ccdb_plasmid'
    elif ccdb_coordinate != (-1,-1) and n20_coordinate ==(-1,-1):
        check_plasmid_label(gb_path, selected_feature_type='misc_feature', target_gene_label=ccdb_label)
        return 'no_sgRNA_plasmid'
    else:
        raise ValueError('The plasmid you uploaded does not contain the necessary tags')

def check_enzyme(enzyme,enzyme_df):
    name = enzyme['enzyme_name']
    gap_seq = enzyme['gap_sequence']
    protection_seq = enzyme['protection_sequence']
    errorMessage = "" 

    if len(gap_seq) < 1 or len(gap_seq) >10:
        errorMessage = f"The gap sequence:{gap_seq} length should be greater than 1 and less than 10" 
        raise ValueError(errorMessage) 
    if not su.is_dna(gap_seq):
        errorMessage = f'The gap sequence:{gap_seq} is not a DNA sequence'
        raise ValueError(errorMessage) 
    if not su.is_dna(protection_seq):
        errorMessage = f'The protection sequence:{gap_seq} is not a DNA sequence'
        raise ValueError(errorMessage) 
    if len(protection_seq) < 1 or len(protection_seq) >20:
        errorMessage = f"The protection sequence:{gap_seq} length should be greater than 1 and less than 10" 
        raise ValueError(errorMessage) 
    
    if name in list(enzyme_df['name']) and len(enzyme_df[enzyme_df['name']==name]) == 1:
        enzyme_df = enzyme_df[enzyme_df['name']==name]
        enzyme_df['protective_base'] = protection_seq
        enzyme_df['gap_len'] = len(gap_seq)
        enzyme_df['gap_seq'] = gap_seq   
        return enzyme_df


def blastEvaluate_uha_dha_in_genome(genome, uha_dha_df, workdir, id, name1, name2, name):

    # run blast   
    workdir = workdir + '/blast/'
  
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    blast_input_file_path = workdir+'/blast.fasta'
    blast_output_file_path= workdir+'/blast_output.txt'
    
    #convert df  
    df = su.convert_twoColumns_to_oneColumns(uha_dha_df, id=id, name1=name1, name2=name2, name=name)
    
    #create fasta
    su.convert_df_to_fastaFile(df, id, name, blast_input_file_path)  
    
    #run blast
    su.blast_ha(genome,blast_input_file_path, blast_output_file_path)

    su.blast_output_evaluate(workdir, blast_input = blast_input_file_path, blast_output= blast_output_file_path)


def create_plasmid_order(uha_dha_sgRNA_df, name, plasmid_type):
        


        if plasmid_type == 'one_plasmid_system':
            plasmid_order_df =  uha_dha_sgRNA_df[['Name', 'Region', 'UHA', 'DHA', 'seq_altered', 'promoter_seq', 'plasmid', 'Target sequence']]
            plasmid_order_df['ID'] = plasmid_order_df['Name'] + ';' + plasmid_order_df['Region']
            plasmid_order_df['SEQUENCE'] = plasmid_order_df['UHA'] + plasmid_order_df['DHA'] + plasmid_order_df['promoter_seq'] + plasmid_order_df['Target sequence']
            
        elif plasmid_type == 'ccdb':
            plasmid_order_df =  uha_dha_sgRNA_df[['Name', 'Region', 'UHA', 'DHA', 'seq_altered', 'plasmid']]
            plasmid_order_df['ID'] = plasmid_order_df['Name'] + ';' + plasmid_order_df['Region']
            plasmid_order_df['SEQUENCE'] = plasmid_order_df['UHA'] + plasmid_order_df['DHA']
           
        elif plasmid_type == 'sgRNA':
            plasmid_order_df =  uha_dha_sgRNA_df[['Name', 'Region',  'promoter_seq', 'Target sequence', 'plasmid']]
            plasmid_order_df['ID'] = plasmid_order_df['Name'] + ';' + plasmid_order_df['Region']
            plasmid_order_df['SEQUENCE'] = plasmid_order_df['promoter_seq'] + plasmid_order_df['Target sequence']

        plasmid_order_df['SEQUENCE SIZE(bp)'] = plasmid_order_df.SEQUENCE.apply(lambda x: len(x))
        plasmid_order_df['IS OR NOT OPTIMIZED'] = 'NOT'  
        plasmid_order_df['PLASMID VECTOR'] = name
        plasmid_order_df['PLASMID SEQUENCE'] = plasmid_order_df['plasmid']
        plasmid_order_df['PLASMID SEQUENCE SIZE(bp)'] = plasmid_order_df['PLASMID SEQUENCE'].apply(lambda x: len(x))
        plasmid_order_df = plasmid_order_df[['ID','SEQUENCE','SEQUENCE SIZE(bp)','IS OR NOT OPTIMIZED', 'PLASMID VECTOR', 'PLASMID SEQUENCE','PLASMID SEQUENCE SIZE(bp)']]
            
        return plasmid_order_df


def add_WT_genome_sequencing_primer_df(uha_dha_sgRNA_df, genome_sequencing_primer_df):

    temp =  uha_dha_sgRNA_df[['Name', 'Region', 'UHA', 'DHA', 'seq_altered', 'ref']]
    temp['ID'] = temp['Name'] +';'+ temp['Region']
    temp = temp[['ID', 'UHA', 'DHA', 'seq_altered', 'ref']]
    columns_list = list(genome_sequencing_primer_df.columns)
    genome_sequencing_primer_df = pd.merge(temp,genome_sequencing_primer_df)
    genome_sequencing_primer_df.fillna('',inplace=True)
    def work(SEQUENCING_TARGET,UHA,ref):

        if SEQUENCING_TARGET == '':
            
            start = SEQUENCING_TARGET.find(UHA)

            before_target = SEQUENCING_TARGET[:start + len(UHA)]
            after_target = SEQUENCING_TARGET[start + len(UHA):]
            
            if ref != '-':
                WT = before_target + ref + after_target
            else:
                WT = SEQUENCING_TARGET
        else:
            WT = SEQUENCING_TARGET
        return WT
    
    genome_sequencing_primer_df['WT_SEQUENCING_TARGET'] = genome_sequencing_primer_df.apply(lambda x: work(x['SEQUENCING_TARGET'], x['UHA'], x['ref']), axis=1)
    genome_sequencing_primer_df['WT_SEQUENCING_TARGET_SIZE'] = genome_sequencing_primer_df.WT_SEQUENCING_TARGET.apply(lambda x: len(x))

    columns_list.extend(['WT_SEQUENCING_TARGET','WT_SEQUENCING_TARGET_SIZE'])
    genome_sequencing_primer_df = genome_sequencing_primer_df[columns_list]

    return genome_sequencing_primer_df

  
def clean_noUsefulContent(path_1, noUsefulContent ='BASE COUNT'):
    
    with open(path_1,"r",encoding="utf-8") as f:
        lines = f.readlines()

    with open(path_1,"w",encoding="utf-8") as f_w:
        for line in lines:
            if noUsefulContent in line:
                continue
            f_w.write(line)

    print(path_1,'过滤完毕！')


def main(data):
      
    chopchop_input = data['chopchop_input']
    print(os.getcwd())
    
    #uha_dha参数
    uha_dha_params = data['uha_dha_config']
    config.UHA_DHA_CONFIG = uha_dha_params

    
    # enzyme_path = parent_base_path +'/'+ data['enzyme_path']
    sgRNA_result_path = data['sgRNA_result_path']
    plasmid_file_1 = data['one_plasmid_file_path']
    plasmid_file_2 = data['no_ccdb_plasmid']
    plasmid_file_3 = data['no_sgRNA_plasmid']
    
    #plasmid label
    ccdb_label = data['plasmid_label']['ccdb_label']
    promoter_terminator_label = data['plasmid_label']['promoter_terminator_label']
    n_20_label = data['plasmid_label']['n_20_label']
    promoter_label = data['plasmid_label']['promoter_label']

    #检查质粒标签
    if ccdb_label == promoter_terminator_label or ccdb_label == n_20_label:
        errorMessage = 'There are duplicates in the three plasmid label you entered'
        raise ValueError(errorMessage)
    elif promoter_terminator_label == n_20_label or promoter_terminator_label == ccdb_label:
        errorMessage = 'There are duplicates in the three plasmid label you entered'
        raise ValueError(errorMessage)
  
   
    one_plasmid_file_path=''
    no_ccdb_plasmid=''
    no_sgRNA_plasmid=''
    
    #0.检查质粒,自动判断质粒类型
    if plasmid_file_1 != '':
        clean_noUsefulContent(plasmid_file_1, noUsefulContent ='BASE COUNT')
        plasmid_type = check_plasmid(plasmid_file_1, ccdb_label, promoter_terminator_label, n_20_label)
        if plasmid_type == 'one_plasmid_file_path':
            one_plasmid_file_path = plasmid_file_1
        elif plasmid_type == 'no_ccdb_plasmid':
            no_ccdb_plasmid = plasmid_file_1
        elif plasmid_type == 'no_sgRNA_plasmid':
            no_sgRNA_plasmid = plasmid_file_1
        elif plasmid_type == 'error':
            errorMessage = 'There is a problem with the plasmid you uploaded'
            raise ValueError(errorMessage)
    
    if plasmid_file_2 != '':
        clean_noUsefulContent(plasmid_file_2, noUsefulContent ='BASE COUNT')
        plasmid_type = check_plasmid(plasmid_file_2, ccdb_label, promoter_terminator_label, n_20_label)
        if plasmid_type == 'one_plasmid_file_path':
            one_plasmid_file_path = plasmid_file_2
        elif plasmid_type == 'no_ccdb_plasmid':
            no_ccdb_plasmid = plasmid_file_2
        elif plasmid_type == 'no_sgRNA_plasmid':
            no_sgRNA_plasmid = plasmid_file_2
        elif plasmid_type == 'error':
            errorMessage = 'There is a problem with the plasmid you uploaded'
            raise ValueError(errorMessage)

    if plasmid_file_3 != '':  
        clean_noUsefulContent(plasmid_file_3, noUsefulContent ='BASE COUNT')
        plasmid_type = check_plasmid(plasmid_file_3, ccdb_label, promoter_terminator_label, n_20_label)    
        if plasmid_type == 'one_plasmid_file_path':
            one_plasmid_file_path = plasmid_file_3
        elif plasmid_type == 'no_ccdb_plasmid':
            no_ccdb_plasmid = plasmid_file_3
        elif plasmid_type == 'no_sgRNA_plasmid':
            no_sgRNA_plasmid = plasmid_file_3
        elif plasmid_type == 'error':
            errorMessage = 'There is a problem with the plasmid you uploaded'
            raise ValueError(errorMessage)


    #1.获取引物，检查引物
    sgRNA_primer_json = data['sgRNA_primer_json']
    ccdb_primer_json = data['ccdb_primer_json']  
    primer_json = data['primer_json']  

    if primer_json != {} and one_plasmid_file_path !='':
        #取质粒骨架，质粒label：ccdb_label、promoter_terminator_label、n_20_label
        plasmid_backbone =  p_d_seq.get_plasmid_backbone_by_labels(one_plasmid_file_path, ccdb_label=ccdb_label, promoter_terminator_label=promoter_terminator_label, n_20_label=n_20_label)
        failture_seq_json = check_quality_control(plasmid_backbone, primer_json)
        if failture_seq_json != {} or failture_seq_json == None:
            raise ValueError(failture_seq_json)
    
    if sgRNA_primer_json !={} and no_ccdb_plasmid !='':
        plasmid_backbone =  p_d_seq.get_plasmid_backbone_by_labels(no_ccdb_plasmid, ccdb_label=ccdb_label, promoter_terminator_label=promoter_terminator_label, n_20_label=n_20_label)
        failture_seq_json = check_quality_control(plasmid_backbone, sgRNA_primer_json)
        if failture_seq_json != {}:
            raise ValueError(failture_seq_json)
            
    if ccdb_primer_json !={} and no_sgRNA_plasmid != '':
        plasmid_backbone =  p_d_seq.get_plasmid_backbone_by_labels(no_sgRNA_plasmid, ccdb_label=ccdb_label, promoter_terminator_label=promoter_terminator_label, n_20_label=n_20_label)
        failture_seq_json = check_quality_control(plasmid_backbone, ccdb_primer_json)
        if failture_seq_json != {}:
            raise ValueError(failture_seq_json)


    #2.根据获取的label，取出region序列，且检查区域序列
    region_label = data['region_label']  
    sgRNA_region_label = data['sgRNA_region_label']
    ccdb_region_label = data['ccdb_region_label'] 

    region_seq_json = {}
    sgRNA_region_seq_json = {}
    ccdb_region_seq_json = {}
    
    if region_label !="" and one_plasmid_file_path !="":
        seq = su.get_sequence_by_feature_label(one_plasmid_file_path, region_label)
        if seq == None:
            raise ValueError('This plasmid feature cannot be sequenced on the plasmid')
        else:
            region_seq_json = {'region1':str(seq)}
            #检查
            plasmid_backbone =  p_d_seq.get_plasmid_backbone_by_labels(one_plasmid_file_path, ccdb_label=ccdb_label, promoter_terminator_label=promoter_terminator_label, n_20_label=n_20_label)
            failture_seq_json = check_quality_control(plasmid_backbone, region_seq_json)
            if failture_seq_json != {}:
                raise ValueError(failture_seq_json)
    
    if sgRNA_region_label !="" and no_ccdb_plasmid !="":
        seq = su.get_sequence_by_feature_label(no_ccdb_plasmid, sgRNA_region_label)
        if seq == None:
            raise ValueError('This plasmid feature cannot be sequenced on the plasmid')
        else:  
            sgRNA_region_seq_json = {'region1':str(seq)}
            #检查
            plasmid_backbone =  p_d_seq.get_plasmid_backbone_by_labels(no_ccdb_plasmid, ccdb_label=ccdb_label, promoter_terminator_label=promoter_terminator_label, n_20_label=n_20_label)
            failture_seq_json = check_quality_control(plasmid_backbone, sgRNA_region_seq_json)
            if failture_seq_json != {}:
                raise ValueError(failture_seq_json)
    
    if ccdb_region_label !="" and no_sgRNA_plasmid !="":
        seq = su.get_sequence_by_feature_label(no_sgRNA_plasmid, ccdb_region_label)
        if seq == None:
            raise ValueError('This plasmid feature cannot be sequenced on the plasmid')
        else:
            ccdb_region_seq_json = {'region1':str(seq)}
            #检查
            plasmid_backbone =  p_d_seq.get_plasmid_backbone_by_labels(no_sgRNA_plasmid, ccdb_label=ccdb_label, promoter_terminator_label=promoter_terminator_label, n_20_label=n_20_label)
            failture_seq_json = check_quality_control(plasmid_backbone, ccdb_region_seq_json)
            if failture_seq_json != {}:
                raise ValueError(failture_seq_json)

  
    #genome
    genome_path = data['ref_genome'] 

    #配置引物参数
    # config.S_GLOBAL_ARGS = data['S_GLOBAL_ARGS']
    # config.Q_ARGS = data['Q_ARGS']
    config.UHA_ARGS = data['UHA_ARGS']
    config.SEQ_ALTERED_ARGS = data['SEQ_ALTERED_ARGS']
    config.DHA_ARGS = data['DHA_ARGS']
    
    #4.判断使用何种技术方法：PCR、OLIGO
    if data.get('UP_SGRNA_ARGS').get('PRIMER_OPT_TM') == '' or data.get('UP_SGRNA_ARGS') == {} or data.get('DOWN_SGRNA_ARGS') == {}:
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

    
    # 5.read 编辑序列信息,给chopchop输入加uha、dha信息
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
    #6.检查酶的相关信息
    enzyme_df = check_enzyme(enzyme, enzyme_df)   


    # 7.判断不同的场景，提取用户选择的sgRNA     
    if scene == 'only_primer': 
        sgRNA = info_input_df[['name','crrna','region']].rename(columns={'crrna':"Target sequence","name":"Name","region":"Region"})
        sgRNA['Rev Target sequence'] = sgRNA['Target sequence'].apply(lambda x: su.revComp(x)) 

    elif scene == 'both_sgRNA_primer':
        selected_sgRNA_result = data['sgRNA_result']
        sgRNA = p_d_seq.extract_sgRNA_from_chopchop(sgRNA_result_path, selected_sgRNA_result)
    

    # 8.设计源生同源臂引物
    uha_dha_primer_df, failture_uha_dha_primer_df = extract_uha_dha_primer(info_input_df, sgRNA)   
    #对uha_dha_primer_df过滤，筛选出成对的UHA，DHA引物
    def work(x):
        if len(x) == 2:
            return x
    uha_dha_primer_df = uha_dha_primer_df.groupby(by='Name').apply(lambda x: work(x))
    uha_dha_primer_df.index = range(0,len(uha_dha_primer_df))
    
    # info_input_df.to_csv('info_input_df.csv',index=False)

    # 9.提取同源臂
    uha_dha_info_primer_df, uha_dha_df, uha_dha_sgRNA_df, info_df = extract_uha_dha(info_input_df, uha_dha_primer_df, sgRNA)
    # uha_dha_df.to_csv('uha_dha_df.csv',index=False)

    # 10.对同源臂在基因组上进行序列比对  
    blastEvaluate_uha_dha_in_genome(genome_path, uha_dha_df, workdir=output, id='Name',name1='UHA',name2='DHA',name='seq')
    
    # 11.判断质粒的执行类型  
    print('判断质粒的执行类型:',one_plasmid_file_path,no_ccdb_plasmid,no_sgRNA_plasmid)  
    if  one_plasmid_file_path != '' and no_ccdb_plasmid == '' and no_sgRNA_plasmid == '':
        plasmid_system_type = 1
    elif  one_plasmid_file_path =='' and no_ccdb_plasmid != '' and no_sgRNA_plasmid != '':
        plasmid_system_type = 2
    elif  one_plasmid_file_path !='' and no_ccdb_plasmid !='' and no_sgRNA_plasmid !='':
        plasmid_system_type = 0
    else:
        raise ValueError("你选择的双质粒系统，质粒没有上传完整!",one_plasmid_file_path,no_ccdb_plasmid,no_sgRNA_plasmid)
        # return '你选择的双质粒系统，质粒没有上传完整!'
    

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
        import time
        start_time = time.time()
        
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
                                                                output,
                                                                ccdb_label,
                                                                promoter_terminator_label,
                                                                n_20_label,
                                                                promoter_label,
                                                                failture_uha_dha_primer_df
                                                                )

        end_time = time.time()
        # 计算函数的执行时间
        exec_time = end_time - start_time
        print("函数的执行时间为：", exec_time, "秒")
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
                                    promoter_label,
                                    output,
                                    failture_uha_dha_primer_df
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
                                                                output,
                                                                ccdb_label,
                                                                promoter_terminator_label,
                                                                n_20_label,
                                                                promoter_label,
                                                                failture_uha_dha_primer_df
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
                                    promoter_label,
                                    output,
                                    failture_uha_dha_primer_df
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
 
    # data1 = {     
    #     "chopchop_input": "./input/only_primer/info_input.csv",   
    #     "sgRNA_result_path": "",
    #     "edit_sequence_design_workdir":"./output/only_primer",
    #     "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
    #     "one_plasmid_file_path":"./input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
    #     "no_ccdb_plasmid":"./input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
    #     "no_sgRNA_plasmid":"./input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
    #     "scene":"only_primer",  
    #     "uha_dha_config": {
    #         "max_right_arm_seq_length": 300,  
    #         "max_left_arm_seq_length": 300,   
    #         "min_left_arm_seq_length": 290,   
    #         "min_right_arm_seq_length": 290     
    #     },

    #     "plasmid_label":{
    #         "ccdb_label":"ccdB",  
    #         "promoter_terminator_label":"gRNA",
    #         "n_20_label":"N20"
    #     },

    #     "primer_json":{
        
    #     },
    #     "region_label":"",       

    #     "sgRNA_primer_json":{

    #     },

    #     "ccdb_primer_json":{
               
    #     },   
    
    #     "sgRNA_region_label":"psiGA1",
        
    #     "ccdb_region_label":"psiGA1",   
        
    #     "enzyme":{
    #         "enzyme_name":"BsaI",
    #         "gap_sequence":"AA",  
    #         "protection_sequence":"CCA"   
    #     },  
          
    #     "UHA_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     "SEQ_ALTERED_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,
    #         "PRIMER_MAX_TM": 75,  
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     "DHA_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,
    #         "PRIMER_MAX_TM": 75,
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #      "UP_SGRNA_ARGS": {
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     "DOWN_SGRNA_ARGS": {
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },

    #     "PLASMID_Q_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     "GENOME_Q_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,     
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     'sgRNA_result':{
    #         "Cgl0006_1176_G_A_sub":"1",
    #         "Cgl2342_213_GCA_ins":"1",
    #         "Cgl1436_1113_CAA_del":"1",
    #         "Cgl1790_1647_TCC_sub":"1",
    #         "Cgl1386_327_18to15_sub":"1",
    #         "Cgl0591_-1_Ppgk_promoter_ins":"1",
    #         "Cgl0141_cds_del":"1",
    #         "153019_ecoil_ybeL_ins":"1",
    #         "Cgl0851_ecoli_pgi_sub":"1"
    #     }      
    # }


    # data2 = {     
    #     "chopchop_input": "./input/only_primer/info_input.csv",   
    #     "sgRNA_result_path": "",
    #     "edit_sequence_design_workdir":"./output/only_primer",
    #     "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
    #     "one_plasmid_file_path":"./input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
    #     "no_ccdb_plasmid":"./input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
    #     "no_sgRNA_plasmid":"./input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
    #     "scene":"only_primer",
    #     "uha_dha_config": {
    #         "max_right_arm_seq_length": 1050,  
    #         "max_left_arm_seq_length": 1050,   
    #         "min_left_arm_seq_length": 1000,   
    #         "min_right_arm_seq_length": 1000     
    #     },

    #       "plasmid_label":{
    #         "ccdb_label":"ccdB",  
    #         "promoter_terminator_label":"gRNA",
    #         "n_20_label":"N20"
    #     },

    #     "primer_json":{   
          
    #     },
    #     "region_label":"",     

    #     "sgRNA_primer_json":{
            
    #     },
    #     "ccdb_primer_json":{
                
    #     },   
    
    #     "sgRNA_region_label":"",
        
    #     "ccdb_region_label":"",   
        
    #     "enzyme":{
    #         "enzyme_name":"BbsI",
    #         "gap_sequence":"AA",    
    #         "protection_sequence":"CCA"   
    #     },      
        
    #     "UHA_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     "SEQ_ALTERED_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },
    #     "DHA_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,
    #         "PRIMER_MAX_TM": 75,
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },

    #     "PLASMID_Q_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },    
    #     "GENOME_Q_ARGS":{
    #         "PRIMER_OPT_TM": 65,
    #         "PRIMER_MIN_TM": 55,     
    #         "PRIMER_MAX_TM": 75,    
    #         "PRIMER_MIN_GC": 20,
    #         "PRIMER_MAX_GC": 80
    #     },  
    #     "UP_SGRNA_ARGS": {
    #         "PRIMER_OPT_TM": "",
    #         "PRIMER_MIN_TM": "",  
    #         "PRIMER_MAX_TM": "",    
    #         "PRIMER_MIN_GC": "",
    #         "PRIMER_MAX_GC": ""
    #     },
    #     "DOWN_SGRNA_ARGS": {
    #         "PRIMER_OPT_TM": "",
    #         "PRIMER_MIN_TM": "",  
    #         "PRIMER_MAX_TM": "",    
    #         "PRIMER_MIN_GC": "",
    #         "PRIMER_MAX_GC": ""
    #     },
    #     'sgRNA_result':{
    #         # "Cgl0006_1176_G_A_sub":"1",
    #         # "Cgl2342_213_GCA_ins":"1",
    #         # "Cgl1436_1113_CAA_del":"1",
    #         # "Cgl1790_1647_TCC_sub":"1",
    #         # "Cgl1386_327_18to15_sub":"1",
    #         # "Cgl0591_-1_Ppgk_promoter_ins":"1",
    #         # "Cgl0141_cds_del":"1",
    #         # "153019_ecoil_ybeL_ins":"1",
    #         # "Cgl0851_ecoli_pgi_sub":"1"
    #     }      
    # }

    # data3= {     
    #     "chopchop_input": "/home/yanghe/tmp/data_preprocessing/output/info_input.csv",   
    #     "sgRNA_result_path": "/home/yanghe/tmp/chopchop/output/sgRNA.csv",
    #     "edit_sequence_design_workdir":"./output/only_primer",
    #     "ref_genome": "/home/yanghe/tmp/data_preprocessing/output/参考基因组eco.fna",
    #     # "one_plasmid_file_path":"/home/yanghe/program/edit_sequence_design/input/only_primer/pMB1-Red-sgRNA-sacB-2297.gb",   
    #     # "one_plasmid_file_path":"/home/yanghe/program/edit_sequence_design/input/Plasmid.gb",
    #     "one_plasmid_file_path":"/home/yanghe/program/edit_sequence_design/input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
    #     "no_ccdb_plasmid":"",
    #     "no_sgRNA_plasmid":"",
    #     "scene":"both_sgRNA_primer",    
    #     "uha_dha_config": {
    #         "max_right_arm_seq_length": 145,  
    #         "max_left_arm_seq_length": 145,   
    #         "min_left_arm_seq_length": 145,   
    #         "min_right_arm_seq_length": 145     
    #     },   

    #     "plasmid_label":{
    #         "ccdb_label":"ccdB",  
    #         "promoter_terminator_label":"gRNA",
    #         "n_20_label":"N20"
    #     },

    #     "primer_json":{
        
    #     },
    #     "region_label":"",       

    #     "sgRNA_primer_json":{
            
    #     },
    #     "ccdb_primer_json":{
            
    #     },   
    
    #     "sgRNA_region_label":"",
            
    #     "ccdb_region_label":"",   
        
    #     "enzyme":{
    #         "enzyme_name":"BsaI",
    #         "gap_sequence":"A",  
    #         "protection_sequence":"CCA"   
    #     },  
          
    #     "UHA_ARGS":{
    #         "PRIMER_OPT_TM": 60,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 65,    
    #         "PRIMER_MIN_GC": 30,
    #         "PRIMER_MAX_GC": 70,
    #         'PRIMER_MIN_SIZE':15,
    #         'PRIMER_MAX_SIZE':25,
    #         # 'PRIMER_OPT_SIZE':20,   
    #     },   
    #     "SEQ_ALTERED_ARGS":{
    #         "PRIMER_OPT_TM": 60,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 65,    
    #         "PRIMER_MIN_GC": 30,
    #         "PRIMER_MAX_GC": 70,
    #         'PRIMER_MIN_SIZE':15,
    #         'PRIMER_MAX_SIZE':25,
    #         # 'PRIMER_OPT_SIZE':20,   
    #     },
    #     "DHA_ARGS":{
    #         "PRIMER_OPT_TM": 60,
    #         "PRIMER_MIN_TM": 55,  
    #         "PRIMER_MAX_TM": 65,    
    #         "PRIMER_MIN_GC": 30,
    #         "PRIMER_MAX_GC": 70,
    #         'PRIMER_MIN_SIZE':15,
    #         'PRIMER_MAX_SIZE':25,
    #         # 'PRIMER_OPT_SIZE':20,   
    #     },
    #      "UP_SGRNA_ARGS": {
    #       "PRIMER_OPT_TM": 60,
    #         "PRIMER_MIN_TM": 40,  
    #         "PRIMER_MAX_TM": 70,    
    #         "PRIMER_MIN_GC": 30,
    #         "PRIMER_MAX_GC": 70,
    #         'PRIMER_MIN_SIZE':14,
    #         'PRIMER_MAX_SIZE':23,
    #         # 'PRIMER_OPT_SIZE':18 
    #     },
    #     "DOWN_SGRNA_ARGS": {
    #         "PRIMER_OPT_TM": 60,
    #         "PRIMER_MIN_TM": 40,  
    #         "PRIMER_MAX_TM": 70,    
    #         "PRIMER_MIN_GC": 30,
    #         "PRIMER_MAX_GC": 70,
    #         'PRIMER_MIN_SIZE':14,
    #         'PRIMER_MAX_SIZE':23,
    #         # 'PRIMER_OPT_SIZE':18,   
    #     },

    #     "PLASMID_Q_ARGS":{
    #         "PRIMER_OPT_TM": 10,
    #         "PRIMER_MIN_TM": 10,  
    #         "PRIMER_MAX_TM": 10,    
    #         "PRIMER_MIN_GC": 10,
    #         "PRIMER_MAX_GC": 10,
    #         'PRIMER_MIN_SIZE':15,
    #         'PRIMER_MAX_SIZE':25,
    #         'PRIMER_OPT_SIZE':20,   
    #     },
    #     "GENOME_Q_ARGS":{
    #         "PRIMER_OPT_TM": 10,
    #         "PRIMER_MIN_TM": 10,  
    #         "PRIMER_MAX_TM": 10,    
    #         "PRIMER_MIN_GC": 10,
    #         "PRIMER_MAX_GC": 10,
    #         'PRIMER_MIN_SIZE':15,
    #         'PRIMER_MAX_SIZE':25,
    #         # 'PRIMER_OPT_SIZE':20,   
    #     },
    #     'sgRNA_result':{
    #         'b4413':"1",
    #         'b4762':"1",
    #         'b4577':"1",
    #         'b4810':"1",
    #         'b4414':"1",
    #         'b4809':"1",
    #         'b4690':"1",
    #         'b0455':"1",
    #         'b4585':"1",
    #         'b4831':"1",
    #         'b4835':"1",
    #         'b4763':"1",
    #         'b4808':"1",
    #         'b4764':"1",
    #         'b4814':"1",
    #         'b4416':"1",
    #         'b4417':"1",
    #         'b4826':"1",
    #         'b4418':"1",
    #         'b4832':"1",
    #         'b4806':"1",
    #         'b4420':"1",
    #         'b4422':"1",
    #         'b4424':"1",
    #         'b4813':"1",
    #         'b4425':"1",
    #         'b4426':"1",
    #         'b4699':"1",
    #         'b4714':"1",
    #         'b4833':"1",
    #         'b4427':"1",
    #         'b4597':"1",
    #         'b4429':"1",
    #         'b4698':"1",
    #         'b1574':"1",
    #         'b4430':"1",
    #         'b4431':"1",
    #         'b4432':"1",
    #         'b4433':"1",
    #         'b4717':"1",
    #         'b4828':"1",
    #         'b4759':"1",
    #         'b4719':"1",
    #         'b4827':"1",
    #         'b1954':"1",
    #         'b4603':"1",
    #         'b4435':"1",
    #         'b4436':"1",
    #         'b4437':"1",
    #         'b4438':"1",
    #         'b4439':"1",
    #         'b4761':"1",
    #         'b4811':"1",
    #         'b4440':"1",
    #         'b4441':"1",
    #         'b4608':"1",
    #         'b4609':"1",
    #         'b4805':"1",
    #         'b2621':"1",
    #         'b4442':"1",
    #         'b4701':"1",
    #         'b4408':"1",
    #         'b4443':"1",
    #         'b4444':"1",
    #         'b4445':"1",
    #         'b2911':"1",
    #         'b4446':"1",
    #         'b4447':"1",
    #         'b4611':"1",
    #         'b3123':"1",
    #         'b4449':"1",
    #         'b4450':"1",
    #         'b4451':"1",
    #         'b4712':"1",
    #         'b4713':"1",
    #         'b4704':"1",
    #         'b4718':"1",
    #         'b4452':"1",
    #         'b4454':"1",
    #         'b4834':"1",
    #         'b4760':"1",
    #         'b4829':"1",
    #         'b4616':"1",
    #         'b4804':"1",
    #         'b4456':"1",
    #         'b4707':"1",
    #         'b3864':"1",
    #         'b4457':"1",
    #         'b4716':"1",
    #         'b4458':"1",
    #         'b4691':"1",
    #         'b4830':"1",
    #         'b4807':"1",
    #         'b4758':"1",
    #         'b4459':"1",
    #         'b4825':"1",
    #         'b4624':"1",
    #         'b4625':"1"
    #     }      
    # }
    
    data1 = {     
        "chopchop_input": "/home/yanghe/tmp/data_preprocessing/output/info_input.csv",   
        "sgRNA_result_path": "",
        "edit_sequence_design_workdir":"/home/yanghe/tmp/edit_sequence_design/output/",
        "ref_genome":"/home/yanghe/tmp/data_preprocessing/output/eco.fna",

        "one_plasmid_file_path":"./input/only_primer/大肠图谱-gRNA反向1-.gb",  
        "no_ccdb_plasmid":"",  
        "no_sgRNA_plasmid":"",

        "scene":"only_primer",

        'sgRNA_result':{}, 

        "uha_dha_config": {
            "max_right_arm_seq_length": 145,  
            "max_left_arm_seq_length": 145,     
            "min_left_arm_seq_length": 145,   
            "min_right_arm_seq_length": 145     
        },
        "plasmid_label":{
            "ccdb_label":"HR",  
            "promoter_terminator_label":"gRNA",
            "n_20_label":"N20",
            "promoter_label":"promoter"
        },
        
        "primer_json":{},
        "region_label":"", 
   
        "sgRNA_primer_json":{},
        "ccdb_primer_json":{},   
        "sgRNA_region_label":"",
        "ccdb_region_label":"",  

        
        "enzyme":{
            "enzyme_name":"BbsI",
            "gap_sequence":"AA",    
            "protection_sequence":"CCA"   
        },      
        
        "UHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },
        "SEQ_ALTERED_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },
        "DHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },

        "PLASMID_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },    
        "GENOME_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,     
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },  
        "UP_SGRNA_ARGS": {
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },
        "DOWN_SGRNA_ARGS": {
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18, 
        },
    }

    data2 = {     
        "chopchop_input": "/home/yanghe/tmp/data_preprocessing/output/info_input.csv",   
        "sgRNA_result_path": "/home/yanghe/tmp/chopchop/output/sgRNA.csv",
        "edit_sequence_design_workdir":"/home/yanghe/tmp/edit_sequence_design/output/",
        "ref_genome":"/home/yanghe/tmp/data_preprocessing/output/xxx.fna",

        "one_plasmid_file_path":"./input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
        "no_ccdb_plasmid":"./input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "no_sgRNA_plasmid":"./input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",

        "scene":"both_sgRNA_primer", 

        "uha_dha_config": {
            "max_right_arm_seq_length": 300,  
            "max_left_arm_seq_length": 300,   
            "min_left_arm_seq_length": 290,   
            "min_right_arm_seq_length": 290     
        },

        "plasmid_label":{
            "ccdb_label":"ccdB",  
            "promoter_terminator_label":"gRNA",
            "n_20_label":"N20",
            "promoter_label":"promoter"
        },

        "primer_json":{
        
        },
        "region_label":"",       

        "sgRNA_primer_json":{

        },

        "ccdb_primer_json":{
               
        },   
    
        "sgRNA_region_label":"",
        
        "ccdb_region_label":"",   
        
        "enzyme":{
            "enzyme_name":"BsaI",
            "gap_sequence":"AA",  
            "protection_sequence":"CCA"   
        },  
          
        "UHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18
        },
        "SEQ_ALTERED_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,  
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18
        },
        "DHA_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,
            "PRIMER_MAX_TM": 75,
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18 
        },
         "UP_SGRNA_ARGS": {
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18
        },
        "DOWN_SGRNA_ARGS": {
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18
        },

        "PLASMID_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18
        },
        "GENOME_Q_ARGS":{
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,     
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            'PRIMER_OPT_GC':65,
            "PRIMER_MAX_GC": 80,
            'PRIMER_MIN_SIZE':15,
            'PRIMER_MAX_SIZE':25,
            'PRIMER_OPT_SIZE':18
            
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

    a=main(data1)     

    print(a) 