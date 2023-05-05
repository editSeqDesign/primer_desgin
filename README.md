# editing sequence design
## Project Introduction 
This project is used for primer design and plasmid map generation to complete biological experiments based on Crispr-HR genetic modification.This system supports biological experiments with single plasmid and double plasmid systems, respectively.

    The application scenario is divided into i.only_primer，ii. both_sgRNA_primer

    i.only_primer
    data1 = {     
        "chopchop_input": ".input/info_input.csv",   
        "sgRNA_result_path": "",
        "edit_sequence_design_workdir":"/home/XXX/tmp/edit_sequence_design/output/",
        "ref_genome":"./input/xxx.fna",
        "one_plasmid_file_path":"./input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
        "no_ccdb_plasmid":"./input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "no_sgRNA_plasmid":"./input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
        "scene":"only_primer",

        'sgRNA_result':{},  
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
            "PRIMER_OPT_TM": "",
            "PRIMER_MIN_TM": "",  
            "PRIMER_MAX_TM": "",    
            "PRIMER_MIN_GC": "",
            "PRIMER_MAX_GC": ""
        },
        "DOWN_SGRNA_ARGS": {
            "PRIMER_OPT_TM": "",
            "PRIMER_MIN_TM": "",  
            "PRIMER_MAX_TM": "",    
            "PRIMER_MIN_GC": "",
            "PRIMER_MAX_GC": ""
        },      
    }
    output:('/home/yanghe/tmp/edit_sequence_design/output/one_plasmid_system_result.zip', '/home/yanghe/tmp/edit_sequence_design/output/two_plasmid_system_result.zip')


    ii.both_sgRNA_primer
    data2 = {     
        "chopchop_input": "/home/XXX/tmp/data_preprocessing/output/info_input.csv",   
        "sgRNA_result_path": "/home/XXX/tmp/chopchop/output/sgRNA.csv",
        "edit_sequence_design_workdir":""/home/XXX/tmp/edit_sequence_design/output/"",
        "ref_genome":"/home/XXX/tmp/data_preprocessing/output/xxx.fna",

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
            "n_20_label":"N20"
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
         "UP_SGRNA_ARGS": {
            "PRIMER_OPT_TM": 65,
            "PRIMER_MIN_TM": 55,  
            "PRIMER_MAX_TM": 75,    
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80
        },
        "DOWN_SGRNA_ARGS": {
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
    output:('/home/yanghe/tmp/edit_sequence_design/output/one_plasmid_system_result.zip', '/home/yanghe/tmp/edit_sequence_design/output/two_plasmid_system_result.zip')

## Installation

```shell
# python3.8

pip install -r requirements.txt

# blastn makeblastdb 需要安装并添加到环境变量
```

## Usage

```shell

git clone {}
cd edit_sequence_design

python main.py

```