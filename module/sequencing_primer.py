import sys,os
# sys.path.append('../')  
from sgRNA_utils.sgRNA_primer_config import config   
import sgRNA_utils.sgRNA_primer_util as util
import pandas as pd
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
from Bio import SeqIO



#生成测序引物模板
def design_sequencing_primers(plasmid_seq,path,dict_plasmid_id, dict_plasmid_seq,mute_position=0,seq_type=''):
    '''
        生成测序引物模板（在target_gene_seq 前后加200bp）,设计测序引物
        params:     
            dict_plasmid_id  
            dict_plasmid_seq
        returns:
    '''
    #生成引物模板
    target_gene_down_seq = dict_plasmid_seq['target_gene_down_seq']
    target_gene_up_seq = dict_plasmid_seq['target_gene_up_seq']
    target_gene_seq = dict_plasmid_seq['mute_after_target_gene_seq']

    
    # target_gene_seq 前后加200bp    
    if len(target_gene_down_seq) >= 200:
        sequencing_peimers_template = target_gene_seq + target_gene_down_seq[:200]
    else:
        sequencing_peimers_template = target_gene_seq + target_gene_down_seq
        temp_len = 200 - len(target_gene_down_seq) 
        sequencing_peimers_template = sequencing_peimers_template + target_gene_up_seq[:temp_len]

    if len(target_gene_up_seq) >= 200:
        sequencing_peimers_template = target_gene_up_seq[-200:] + sequencing_peimers_template
    else:
        sequencing_peimers_template = target_gene_up_seq + sequencing_peimers_template  
        temp_len = 200 - len(target_gene_up_seq)
        sequencing_peimers_template = target_gene_down_seq[-temp_len:] + sequencing_peimers_template
    #设计引物 
    if mute_position ==0:
        result, failture_result = design_seq_primer(plasmid_seq, path,seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq,seq_type=seq_type)
        print(result, failture_result)
    else:
  
         result, failture_result = design_seq_primer(plasmid_seq, path,seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq,mute_position=mute_position,seq_type=seq_type)
         
    return result, failture_result

#设计测序引物
def design_seq_primer(plasmid_seq,path,seq_id, seq_temp, seq_target, mute_position=0,seq_type=''):
    dict_seq_primer={}
    dict_seq_primer_failtrue={}
    len_target=len(seq_target)     
    site_target_temp=seq_temp.find(seq_target)
    temp_p1=seq_temp[site_target_temp-120:site_target_temp-80]



    #第一条引物
    dict_res_p1 = util.primer_design(seqId=seq_id, 
                                        seqTemplate=temp_p1,
                                        stype='none',
                                        mute_type='sequencing',
                                        primer_type=seq_type,
                                        global_args=config.Q_GLOBAL_ARGS
                                    )
    u = judge_primer_is_or_not( plasmid_seq,
                            path,
                            seq_id,
                            dict_res_p1,
                            dict_seq_primer,
                            dict_seq_primer_failtrue,
                            primer_name="SEQUENCING_PRIMER_1",
                            type='LEFT')
 
    if dict_res_p1.get(f'PRIMER_LEFT_{u}') != None and u != -1: 

        if dict_res_p1.get(f'PRIMER_LEFT_{u}') != None :
            SEQUENCING_TARGET_START = site_target_temp - 120 + dict_res_p1[f'PRIMER_LEFT_{u}'][0]
        else:
            SEQUENCING_TARGET_START = -1
    else:
        SEQUENCING_TARGET_START = -1
    
    # if 600 < len_target <= 1200:
    if seq_type == 'genome_seq':

        #第二条引物
        temp_p2 = seq_temp[site_target_temp + len_target + 80 : site_target_temp + len_target + 120]
        dict_res_p2 = util.primer_design(seqId=seq_id,
                                            seqTemplate=temp_p2,
                                            stype='none',
                                            mute_type='sequencing',
                                            primer_type=seq_type,
                                            global_args=config.Q_GLOBAL_ARGS
                                        )
            #判断引物是否成功
        u = judge_primer_is_or_not( plasmid_seq,path,seq_id,dict_res_p2,
                                dict_seq_primer,
                                dict_seq_primer_failtrue,
                                primer_name="SEQUENCING_PRIMER_2",
                                type='RIGHT')
        
        if dict_res_p2.get(f'PRIMER_LEFT_{u}') != None and dict_res_p1.get(f'PRIMER_LEFT_{u}') != None and u != -1:

            if SEQUENCING_TARGET_START != -1:
                print( dict_res_p2[f'PRIMER_LEFT_{u}'][0], dict_res_p2[f'PRIMER_LEFT_{u}'][1])
                SEQUENCING_TARGET_END = site_target_temp + len_target + 80 + dict_res_p2[f'PRIMER_LEFT_{u}'][0] + dict_res_p2[f'PRIMER_LEFT_{u}'][1]
                dict_seq_primer['SEQUENCING_TARGET'] = seq_temp[SEQUENCING_TARGET_START:SEQUENCING_TARGET_END]
                dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = SEQUENCING_TARGET_END - SEQUENCING_TARGET_START
            
        
    elif seq_type=='plasmid_seq':

        if 600 < len_target <= 1200:
            temp_p2 = seq_temp[site_target_temp + len_target + 80 : site_target_temp + len_target + 120]
            dict_res_p2=util.primer_design(seqId=seq_id,
                                                seqTemplate=temp_p2,
                                                stype='none',
                                                mute_type='sequencing',
                                                primer_type=seq_type,
                                                global_args=config.Q_GLOBAL_ARGS
                                            )
                #判断引物是否成功
            u = judge_primer_is_or_not( plasmid_seq,path,seq_id,dict_res_p2,
                                    dict_seq_primer,
                                    dict_seq_primer_failtrue,
                                    primer_name="SEQUENCING_PRIMER_2",
                                    type='RIGHT')

    if 1200 < len_target:    
 
        temp_p2 = seq_temp[site_target_temp+len_target+80:site_target_temp+len_target+120]
        dict_res_p2=util.primer_design(seqId=seq_id,
                                        seqTemplate=temp_p2,
                                        stype='none',
                                        mute_type='sequencing',
                                        primer_type=seq_type,
                                        global_args=config.Q_GLOBAL_ARGS
                                    )  
         #判断引物是否成功
        u = judge_primer_is_or_not( plasmid_seq,path,seq_id,dict_res_p2,
                                dict_seq_primer,
                                dict_seq_primer_failtrue,
                                primer_name="SEQUENCING_PRIMER_2",
                                type='RIGHT')

        i=2
        while len_target >= 600:
            temp_pn=seq_temp[site_target_temp+480:site_target_temp+520]
            dict_res_pn = util.primer_design(seqId=seq_id,
                                                seqTemplate=temp_pn,
                                                stype='none',
                                                mute_type='sequencing',
                                                primer_type=seq_type,
                                                global_args=config.Q_GLOBAL_ARGS
                                            )
            seq_target = seq_target[600:]
            site_target_temp=seq_temp.find(seq_target)   
            len_target=len(seq_target)
            if len_target>=600:
                i+=1
                 #判断引物是否成功
                u = judge_primer_is_or_not(plasmid_seq,path,seq_id,
                                        dict_res_pn,
                                        dict_seq_primer,
                                        dict_seq_primer_failtrue,
                                        primer_name=f"SEQUENCING_PRIMER_{i}",
                                        type='LEFT')

        if 200 <= len_target < 600:
            length_add=int((600-len_target)/2)
            temp_pe=seq_temp[site_target_temp-length_add-120:site_target_temp-length_add-80]
            dict_res_pe = util.primer_design(seqId=seq_id,
                                            seqTemplate=temp_pe,
                                            stype='none',
                                            mute_type='sequencing',
                                            primer_type=seq_type,
                                            global_args=config.Q_GLOBAL_ARGS  
                                        )
            #判断引物是否成功
            u = judge_primer_is_or_not( plasmid_seq,path,seq_id,dict_res_pe,
                                    dict_seq_primer,   
                                    dict_seq_primer_failtrue,
                                    primer_name=f"SEQUENCING_PRIMER_{i}",     
                                    type='LEFT')

    return dict_seq_primer, dict_seq_primer_failtrue  
  
#判断引物是与否
def judge_primer_is_or_not(plasmid_seq, path, seq_id, dict_res, primer_suc, primer_fail, primer_name, type='LEFT'):

    #检索引物的个数：
    
   

    if 'b4422' in seq_id:
        print('jhdfsgcjhdgf')

    u = -1
    if path.split('.')[-1] != 'gb':

        plasmid_seq

        record_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))   
        chrom = seq_id.split(';')[1].split(':')[0]
        seq=str(record_dict[chrom].seq)

    else:
        seq = plasmid_seq

    if len(list(dict_res.keys())) < 10:
        primer_fail[f'{seq_id}:{primer_name}'] = dict_res 
    else:

        #引物的数量
        primer_num = dict_res.get('PRIMER_LEFT_NUM_RETURNED')

        #使用blast方法进行序列比对  
        # result_df = blast_primer_in_genome(dict_res, type, path) 
        off_target_dict = {}
        u = 0
        
        for i,k in dict_res.items():
            id = f'PRIMER_{type}_{u}_SEQUENCE'
            primer_seq=dict_res.get(f'PRIMER_{type}_{u}_SEQUENCE')

            if primer_seq != None:  
                j = 0
                #对每一条引物进行脱靶分析，种子序列为15bp
                for m in range(15,len(primer_seq)+1):                     #控制滑动次数
                    if type == 'RIGHT':
                        temp_seq = util.revComp(primer_seq)
                        pattern = temp_seq[j:m]
                    else:
                        pattern = primer_seq[j:m]
                    result = util.find_coordinate_by_pattern(pattern=pattern, seq=seq)
                    j = j + 1

                    #脱靶
                    if result[1] >1:
                        u = u + 1
                        off_target_dict.update({f'{primer_name}_{u}':primer_seq})
                        break
                    #非脱靶
                    elif result[1] == 1:
                        continue
                    else:
                        break
                
                
                if j == len(primer_seq) - 15 + 1  and u < primer_num:      #滑动次数满，且引物数没满，设计成功
                    try:
                        print(seq_id ,dict_res)
                        primer_suc[primer_name] = dict_res[f'PRIMER_{type}_{u}_SEQUENCE']
                        primer_suc[f'{primer_name}_TM'] =dict_res[f'PRIMER_{type}_{u}_TM']
                    except:
                        print(1223333) 
                    break  
                elif  u == primer_num and j < len(primer_seq) - 15 + 1 :        #滑动次数未满，且引物数满，脱靶
                    u = u - 1
                    primer_fail[f'{seq_id}:{primer_name}:off_target'] = off_target_dict
                    break
                    
            else:
               primer_fail[f'{seq_id}:{primer_name}'] = dict_res         #一条引物都没设计出来，设计失败

    return u

#有点问题
def blast_primer_in_genome(dict_res, type, path, ):

    df = pd.DataFrame()
    u = 0 
    for i,k in dict_res.items():
        id = f'PRIMER_{type}_{u}_SEQUENCE'
        seq=dict_res.get(f'PRIMER_{type}_{u}_SEQUENCE')
        if seq != None:
            df = df.append(pd.DataFrame([{'id':id,'seq':seq}]))
            u = u +1

    primer_workdir = config.workdir+'/primer'
    if not exists(primer_workdir):
        os.makedirs(primer_workdir)
        
    blast_output_file_path= os.path.join(primer_workdir,'blast_output.txt')
    blast_input_file_path = os.path.join(primer_workdir,'blast.fasta')
        
    #create fasta
    util.convert_df_to_fastaFile(df,id='id',name='seq',input_fasta = blast_input_file_path)

    #run blast
    if path.split('.')[-1] == 'gb':
        fna_file = config.workdir+'/primer/xxx.fna'
            #convert fasta
        util.gb_2_fna(path, fna_file)
        genome = fna_file
    else:
        genome = path

    util.blast_ha(genome, blast_input_file_path, blast_output_file_path)

    result_df = util.blast_output_evaluate(workdir = config.workdir + '/primer/' , blast_input = blast_input_file_path, blast_output= blast_output_file_path)

    return result_df 









              
  
    

