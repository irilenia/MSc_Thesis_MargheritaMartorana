import configparser
from modules import *

config = configparser.ConfigParser()
config.read('config.ini')

# Importing file names from config.ini file
infile = config.get('Files', 'infile')
temp_file = config.get('Files', 'temp_file')
temp_copy_file = config.get('Files', 'temp_copy_file')
rnafold_fasta_file = config.get('Files', 'rnafold_fasta_file')
rnafold_csv_file = config.get('Files', 'rnafold_csv_file')
rnafold_cluster_fasta_file = config.get('Files', 'rnafold_cluster_fasta_file')
rnafold_cluster_csv_file = config.get('Files', 'rnafold_cluster_csv_file')
final_csv_file = config.get('Files', 'final_csv_file')
final_file = config.get('Files', 'final_file')

# Importing variables from config.ini file
fragment_step = int(config.get('Variables', 'fragment_step'))
fragment_length = int(config.get('Variables', 'fragment_length'))
fragment_overlap = int(config.get('Variables', 'fragment_overlap'))
n_random = int(config.get('Variables', 'n_random'))
pvalue = float(config.get('Variables', 'pvalue'))
energy = float(config.get('Variables', 'energy'))

#Importing commands from config.ini file
cmd_shuffle = config.get('Commands', 'cmd_shuffle')
cmd_fold = config.get('Commands', 'cmd_fold')
cmd_fold_cluster = config.get('Commands', 'cmd_fold_cluster')

def main():
     cmd = ('rm -f -r Output; mkdir Output; mkdir Output/FastaOutput; mkdir Output/CsvOutput')
     process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
     if config.getboolean('Settings', 'transcribe_rna'):
          print('Reading Input File...')
          gene_sequence = doTranscribeRNA(infile)
     if config.getboolean('Settings', 'extract_locations'):
          print('Extracting Gene Location...')
          gene_start_location = getStartLocation(infile)
     if config.getboolean('Settings', 'create_fragments_dictionary'):
          print('Creating Fragments...')
          fragments_dictionary = doFragmentTranscript(gene_sequence, gene_start_location, \
                                 int(fragment_step), int(fragment_length))
     if config.getboolean('Settings', 'shuffle_and_fold'):
          print('Running MEME-shuffle and RNAfold...')
          runShuffleAndFold(infile, fragments_dictionary, temp_file,\
                            temp_copy_file, cmd_shuffle, cmd_fold)
     if config.getboolean('Settings', 'create_csv_file'):
          print('Creating csv file...')
          doCsvFile(rnafold_fasta_file, n_random)
     if config.getboolean('Settings', 'calculate_zscore'):
          print('Calculating Zscore...')
          doZscore(rnafold_csv_file, n_random)
     if config.getboolean('Settings', 'calculate_pvalue'):
          print('Calculating Pvalue...')
          doPvalue(rnafold_csv_file)
     if config.getboolean('Settings', 'create_clusters'):
          print('Creating Clusters...')
          clusters_list = doClusters(rnafold_csv_file, float(pvalue), float(energy), int(fragment_length), \
                                     gene_sequence, gene_start_location, int(fragment_overlap))
     if config.getboolean('Settings', 'shuffle_and_fold_clusters'):
          print('Running MEME-shuffle and RNAfold on clusters...')
          doClusterShuffleAndFold(infile, clusters_list, gene_start_location, gene_sequence, \
                                  rnafold_cluster_fasta_file, temp_file, \
                                  temp_copy_file, cmd_shuffle, cmd_fold_cluster)
     if config.getboolean('Settings', 'create_csv_file_clusters'):
          print('Creating csv file...')
          doCsvFile(rnafold_cluster_fasta_file, n_random)
     if config.getboolean('Settings', 'calculate_zscore_clusters'):
          print('Calculating Zscore...')
          doZscore(rnafold_cluster_csv_file, n_random)
     if config.getboolean('Settings', 'calculate_pvalue_clusters'):
          print('Calculating Pvalue...')
          doPvalue(rnafold_cluster_csv_file)
     if config.getboolean('Settings', 'create_feature_file'):    
          print('Looking for features...')
          doFeatureFile(rnafold_cluster_csv_file, energy, pvalue, \
                        dict_transcripts, final_csv_file, final_file)
     cmd = ('mv *_output.csv Output/CsvOutput; mv *_output.fasta Output/FastaOutput; rm -f temp.fasta')
     process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)    
     print('Program run successfully!')
     

if __name__ == '__main__':
     main()



