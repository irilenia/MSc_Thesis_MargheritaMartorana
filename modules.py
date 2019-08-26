# Import systems
import re, os, sys, csv, scipy, shutil, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats
from Bio.Alphabet import IUPAC
from subprocess import PIPE, run
from transcripts_info import dict_transcripts

def out(command):
     result = run(command, stdout = PIPE, stderr = PIPE, \
                  universal_newlines = True, shell = True)
     return result.stdout

# DCX is a reverse transcribed gene, therefore the start location
# is actually the 'end' location
def getStartLocation(input_file):
     records = SeqIO.parse(open(input_file), 'fasta')
     record = next(records)
     parts = record.description.split(':')
     location = parts[-2]
     if location.isdigit():
          start_location = int(location)
     else:
          start_location = int(len(record.seq))
     return start_location

# Not needed in this case
def getEndLocation(input_file):
     records = SeqIO.parse(open(input_file), 'fasta')
     record = next(records)
     parts = record.description.split(':')
     end_location = parts[-3]
     return end_location

def doTranscribeRNA(input_file):
     for record in SeqIO.parse(open(input_file), 'fasta'):
          sequence = record.seq
     return sequence

def doFragmentTranscript(gene_sequence, gene_start_location, \
                         fragment_step, fragment_length):
     d = {}
     count = 0
     start = gene_start_location
     for i in range(0, len(str(gene_sequence)), fragment_step):
        v = str(gene_sequence)[(count):(count+fragment_length)]
        if len(v) == fragment_length:
             k = str(start) + '_' + str(start-len(v)+1)
             d[k] = v
             count += fragment_step
             start -= fragment_step
        elif len(v) == len(v)%fragment_length:
             break
     return d

def runShuffleAndFold(input_file, dictionary, temp_file,\
                      temp_copy_file, cmd_shuffle, cmd_fold):
     f = input_file.split('.')
     file_name = f[0] + '_output.' + f[1]
     ofile = open(file_name, 'w')
     for key, value in dictionary.items():
          ofile.write('>' + key + '\n' + \
                      value + '\n')
          tempofile = open(temp_file, 'w')
          tempofile.write('>' + key + '\n' + \
                      value + '\n')
          tempofile.close()
          shuffled = out(cmd_shuffle)
          ofile.write(shuffled)
     ofile.close()
     # run RNAfold
     shutil.copyfile(file_name, temp_copy_file)
     in_RNAfold = open(temp_copy_file, 'r')
     out(cmd_fold)
     # remove/move files
     cmd = ('find . -name "*.ps" -delete')
     process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
     os.remove(temp_file)
     os.remove(temp_copy_file)
     return

def doCsvFile(fasta_file, n_random):
     with open (fasta_file, 'r') as f:
          lines = f.read().splitlines()
     # list of IDs
     id_line = lines[0::3]
     secondary_IDs = [] # list of secondary IDs
     for x in id_line:
          x = x.replace('>', '')
          secondary_IDs.append(x)
     item_IDs = secondary_IDs[::int(n_random+1)]
     primary_IDs = [] # list of primary IDs
     for x in item_IDs:
          x = x.replace('>', '')
          count = 0
          while count < (n_random+1):
               primary_IDs.append(x)
               count += 1
     # list of sequences
     sequences = list(lines[1::3])
     # list of fold
     third_line = list(lines[2::3])
     fold = []
     p_1 = re.compile(r'^([\S]+)') # regex for fold
     for i in third_line:
          fold += p_1.findall(i)
     # list of energy
     energy = []
     p_2 = re.compile(r'\s\((.*?)\)') # regex for energy
     for i in third_line:
          energy += p_2.findall(i)
     # make csv from dictionary
     d = {'Primary ID':primary_IDs, 'Secondary ID':secondary_IDs, \
     'Sequence':sequences, 'Fold':fold, 'MFE':energy}
     df = pd.DataFrame(d)
     file_name = (fasta_file.split('.'))[0] + '.csv'
     with open (file_name, 'w') as f:
          df.to_csv(f, sep = '\t', index=False)
          f.close()
     return

def doZscore(csv_file, n_random):
     df = pd.read_csv(csv_file, sep = '\t')
     row_start = 0
     row_end = n_random + 1
     step = n_random + 1
     zscore = []
     while row_end <= len(df):
          selected_rows = df['MFE'].iloc[row_start:row_end]
          arr = []
          for x in selected_rows:
               arr.append(round(float(x)))
          if arr.count(0) == len(arr):
               scores = arr
          else:
               scores = stats.zscore(arr)
          for i in scores:
               zscore.append(round(i, 3))
          arr.clear()
          row_start += step
          row_end += step
     df['Zscore'] = zscore
     with open(csv_file, 'w') as f:
          df.to_csv(f, index=False)
          f.close()
     return


def doPvalue(csv_file):
     df = pd.read_csv(csv_file)
     pvalue = []
     selected_rows = df['Zscore']
     for x in selected_rows:
          pvalue.append(float((round(scipy.stats.norm.sf((abs(x))*2), 3))))
     df['Pvalue'] = pvalue
     with open(csv_file, 'w') as f:
          df.to_csv(f, index=False)
          f.close()
     return

def doClusters(csv_file, pvalue, energy, fragment_length, \
               gene_sequence, gene_start_location, \
               fragment_overlap):
     # Reading first RNA fold output csv file and selecting
     # rows of interest and append them onto a list
     df = pd.read_csv(csv_file)
     df_primary_rows = pd.DataFrame()
     df_primary_rows = df[(df['Primary ID'] == df['Secondary ID']) &\
                          (df['MFE'] <= energy) &\
                          (df['Pvalue'] <= pvalue)]
     primary_list = df_primary_rows['Primary ID'].tolist()
     # Making a list of tuples for the locations
     loc_list = []
     for x in primary_list:
          start_location = (x.split('_'))[1]
          end_location = (x.split('_'))[0]
          loc_list.append((start_location, end_location))
     # Sorting locations
     sorted_by_lower_bound = sorted(loc_list, key=lambda tup: tup[0])
     merged = []
     # Making clusters
     for higher in sorted_by_lower_bound:
          if not merged:
               merged.append(higher)
          else:
               lower = merged[-1]
               if int(lower[1])-int(higher[0]) >= int(fragment_overlap):
                    upper_bound = max(lower[1], higher[1])
                    merged[-1] = (lower[0], upper_bound)
               else:
                    if (int(higher[1])-int(higher[0])) >= int(int(fragment_length)-1):
                         merged.append(higher)
     final_cluster_list = []
     for x in merged:
          if int(x[1])-int(x[0]) > int(int(fragment_length)-1):
               final_cluster_list.append(x)
          else:
               final_cluster_list.append(x)
     final_cluster_list.sort(reverse=True)
     return final_cluster_list 

def doClusterShuffleAndFold(input_file, cluster_list, gene_start_location, \
                            gene_sequence, cluster_fasta_file_name, temp_file, \
                            temp_copy_file, cmd_shuffle, cmd_fold_cluster):
     cluster_dictionary = {}
     for cluster in cluster_list:
          cluster_start = cluster[1]
          cluster_end = cluster[0]
          cluster_location = cluster_end + '_' + cluster_start
          len_cluster = int(cluster_start) - int(cluster_end)
          index_start_sequence = int(gene_start_location) - int(cluster_start)
          index_end_sequence = int(index_start_sequence) + int(len_cluster)
          cluster_sequence = (str(gene_sequence))[slice(index_start_sequence, \
                                                 index_end_sequence)]
          cluster_dictionary[cluster_location] = cluster_sequence
     f = input_file.split('.')
     file_name = f[0] + '_cluster_output.' + f[1]
     ofile = open(file_name, 'w')
     for key, value in cluster_dictionary.items():
          ofile.write('>' + key + '\n' + value + '\n')
          tempofile = open(temp_file, 'w')
          tempofile.write('>' + key + '\n' + value + '\n')
          tempofile.close()
          shuffled = out(cmd_shuffle)
          ofile.write(shuffled)
     ofile.close()
     # Run RNAfold
     shutil.copyfile(file_name, temp_copy_file)
     in_RNAfold = open(temp_copy_file, 'r')
     out(cmd_fold_cluster)
     # remove/move files
     cmd = ('find . -name "*.ps" -delete')
     process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
     os.remove(temp_file)
     os.remove(temp_copy_file)
     return

def doFeatureFile(cluster_csv_file, energy, pvalue, \
                  dict_transcripts, final_csv_file, \
                  final_file):
     # Reading cluster csv file and selecting rows of
     # interest and append them into a list and a list
     # of tuples needed for dicitonary
     df = pd.read_csv(cluster_csv_file)
     df_primary_rows = pd.DataFrame()
     for index, columns in df.iterrows():
          df_primary_rows = df[(df['Primary ID'] == df['Secondary ID']) &\
                               (df['MFE'] <= energy) &\
                               (df['Pvalue'] <= pvalue)]
     cluster_location_list = []
     for index, row in df_primary_rows.iterrows():
          cluster_location_list.append(row['Primary ID'])
     cluster_location_dictionary = {}
     for x in cluster_location_list:
          start_location = (x.split('_'))[1]
          end_location = (x.split('_'))[0]
          cluster_location_dictionary[x] = ((end_location, start_location))
     feature_dict = {}
     count = 0
     while count < 2:
          for loc in cluster_location_dictionary:
               temp_dict = {}
               temp_list = []
               # Iterating through the transcripts dictionary
               for transcript_id, transcript_feature in dict_transcripts.items():
                    # Iterating through the each transcript dictionary
                    for feature in transcript_feature:
                         if int(cluster_location_dictionary[loc][0]) <= int(transcript_feature[feature][0]) and \
                            int(cluster_location_dictionary[loc][0]) >= int(transcript_feature[feature][1]):
                              temp_dict.setdefault(transcript_id, [])
                              if feature not in temp_dict[transcript_id]:
                                   temp_dict[transcript_id].append(feature)
               count += 1
               # Iterating through the transcripts dictionary
               for transcript_id, transcript_feature in dict_transcripts.items():
                    # Iterating through the each transcript dictionary
                    for feature in transcript_feature:
                         if int(cluster_location_dictionary[loc][1]) <= int(transcript_feature[feature][0]) and \
                            int(cluster_location_dictionary[loc][1]) >= int(transcript_feature[feature][1]):
                              temp_dict.setdefault(transcript_id, [])
                              if feature not in temp_dict[transcript_id]:
                                   temp_dict[transcript_id].append(feature)
               for key, value in temp_dict.items():
                    temp_list.append((key, value))
               feature_dict[loc] = temp_list
               count += 1
     # Create final csv file
     df_final = pd.DataFrame()
     df_final['Primary ID'] = df_primary_rows['Primary ID']
     df_final['Sequence'] = df_primary_rows['Sequence']
     df_final['Fold'] = df_primary_rows['Fold']
     df_final['MFE'] = df_primary_rows['MFE']
     df_final['Pvalue'] = df_primary_rows['Pvalue']
     with open (final_csv_file, 'w') as f:
          df_final.to_csv(f, sep='\t', index = False)
          f.close()
     ofile = open(final_file, 'w')
     for key, value in feature_dict.items():
          ofile.write('>' + key + '\n' + str(value) + '\n') 
     return    
