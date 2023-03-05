import os

class config: 
    #单双点设计引物的全局参数  
    GLOBAL_ARGS = {
                'PRIMER_PICK_ANYWAY':1,
                'PRIMER_PRODUCT_SIZE_RANGE': 0,
                'PRIMER_NUM_RETURN':2
        }     

    S_GLOBAL_ARGS = {
            'PRIMER_OPT_SIZE': 20,   
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 65.0,
            'PRIMER_MIN_TM': 55.0,
            'PRIMER_MAX_TM': 75.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
    }
    UHA_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }
    SEQ_ALTERED_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }
    DHA_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }
    UP_SGRNA_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }
    DOWN_SGRNA_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }

    #测序引物设计的全局参数
    Q_ARGS = {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,   
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 65.0,
                'PRIMER_MIN_TM': 55.0,  
                'PRIMER_MAX_TM': 75.0,    
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
    }

    PLASMID_Q_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,  
        'PRIMER_MAX_TM': 75.0,    
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }
    GENOME_Q_ARGS = {
        'PRIMER_OPT_TM': 65.0,
        'PRIMER_MIN_TM': 55.0,  
        'PRIMER_MAX_TM': 75.0,    
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
    }


































    
    Q_GLOBAL_ARGS = {   
                'PRIMER_PICK_ANYWAY':1,    
                'PRIMER_TASK': 'pick_primer_list', 
        }  

    #提取target_gene_seq和marker_seq需要的常量
    AMPLICONIC_GENE_TARGET_SEQ_LENGTH = 0
    
    #提取marker时需要的常量
    AMPLICONIC_MARKER_SEQ_LENGTH = 0

    INPUT_FILE_PATH = ''   
    
    PLASMID_FILE_NAME = ''    
    INPUT_FILE_NAME = ''  
    NEW_OUTPUT_FILE_PATH=''
    OUTPUT_ORDERS_NAME = 'order.xlsx'  
    INPUT_IMG_PATH = '/input/'
    IMG_ORDERS_NAME = 'tib.png'

    SEQUENCING_PRIMER_SUCCESS = 'sequencing_primer_success.xlsx'
    SEQUENCING_PRIMER_FAILTRUE = 'sequencing_primer_failtrue.xlsx'

    OUTPUT_FILE_NAME_PLASMID_MUTATION = "plasmid_mutation.gb" 
    DATA_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) + '/data'

    OUTPUT_FILE_PATH = DATA_ROOT + '/output/'
