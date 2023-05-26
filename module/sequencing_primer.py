import sys,os
# sys.path.append('../')  
from sgRNA_utils.sgRNA_primer_config import config 
import sgRNA_utils.sgRNA_primer_util as util

#生成测序引物模板
def design_sequencing_primers(dict_plasmid_id, dict_plasmid_seq,mute_position=0,seq_type=''):
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
        result, failture_result = design_seq_primer(seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq,seq_type=seq_type)
        print(result, failture_result)
    else:
  
         result, failture_result = design_seq_primer(seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq,mute_position=mute_position,seq_type=seq_type)
         
    return result, failture_result

#设计测序引物
def design_seq_primer(seq_id, seq_temp, seq_target, mute_position=0,seq_type=''):
    dict_seq_primer={}
    dict_seq_primer_failtrue={}
    len_target=len(seq_target)     
    site_target_temp=seq_temp.find(seq_target)
    temp_p1=seq_temp[site_target_temp-120:site_target_temp-80]

    dict_res_p1 = util.primer_design(seqId=seq_id, 
                                        seqTemplate=temp_p1,
                                        stype='none',
                                        mute_type='sequencing',
                                        primer_type=seq_type,
                                        global_args=config.Q_GLOBAL_ARGS
                                    )
   
    #第一条引物
    judge_primer_is_or_not( seq_id,dict_res_p1,   
                            dict_seq_primer,
                            dict_seq_primer_failtrue,
                            primer_name="SEQUENCING_PRIMER_1",
                            type='LEFT')
    
    if dict_res_p1.get('PRIMER_LEFT_0') != None: 
        print(   '--------------------------',dict_res_p1.get('PRIMER_LEFT_0') )
        SEQUENCING_TARGET_START = site_target_temp - 120 + dict_res_p1['PRIMER_LEFT_0'][0]
    
    # if 600 < len_target <= 1200:
    if seq_type == 'genome_seq':

        temp_p2 = seq_temp[site_target_temp + len_target + 80 : site_target_temp + len_target + 120]
        dict_res_p2=util.primer_design(seqId=seq_id,
                                            seqTemplate=temp_p2,
                                            stype='none',
                                            mute_type='sequencing',
                                            primer_type=seq_type,
                                            global_args=config.Q_GLOBAL_ARGS
                                        )
            #判断引物是否成功
        judge_primer_is_or_not( seq_id,dict_res_p2,
                                dict_seq_primer,
                                dict_seq_primer_failtrue,
                                primer_name="SEQUENCING_PRIMER_2",
                                type='RIGHT')
        
        if dict_res_p2.get('PRIMER_LEFT_0') != None and dict_res_p1.get('PRIMER_LEFT_0') != None:
            
            print( dict_res_p2['PRIMER_LEFT_0'][0], dict_res_p2['PRIMER_LEFT_0'][1])
            SEQUENCING_TARGET_END = site_target_temp + len_target + 80 + dict_res_p2['PRIMER_LEFT_0'][0] + dict_res_p2['PRIMER_LEFT_0'][1]

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
            judge_primer_is_or_not( seq_id,dict_res_p2,
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
        judge_primer_is_or_not( seq_id,dict_res_p2,
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
                judge_primer_is_or_not(seq_id,
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
            judge_primer_is_or_not( seq_id,dict_res_pe,
                                    dict_seq_primer,   
                                    dict_seq_primer_failtrue,
                                    primer_name=f"SEQUENCING_PRIMER_{i}",     
                                    type='LEFT')



    return dict_seq_primer, dict_seq_primer_failtrue  
  
#判断引物是与否
def judge_primer_is_or_not( seq_id, dict_res,primer_suc,primer_fail,primer_name,type='LEFT'):




    if len(list(dict_res.keys())) < 10:
        primer_fail[f'{seq_id}:{primer_name}'] = dict_res 
    else:
        if dict_res.get(f'PRIMER_{type}_0_SEQUENCE') == None:
            primer_fail[f'{seq_id}:{primer_name}'] = dict_res 
        else:
            primer_suc[primer_name] = dict_res[f'PRIMER_{type}_0_SEQUENCE']
            primer_suc[f'{primer_name}_TM'] =dict_res[f'PRIMER_{type}_0_TM']
        



              
  
    

