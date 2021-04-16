from os.path import isfile

''' refseq_file_exists
    
    checks if the refseq_summary_file exists in the desired media directory
    
    :returns
        :type boolean
'''
def refseq_file_exists():
    return isfile('media/databases/refseq_summary_file/assembly_summary_refseq.txt')