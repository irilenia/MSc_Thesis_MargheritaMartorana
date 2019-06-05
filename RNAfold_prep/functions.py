import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def getStartLocation(input_file):
     records = SeqIO.parse(open(input_file), 'fasta')
     record = next(records)
     parts = record.description.split(':')
     start_location = int(parts[-3])
     return start_location
     

def doFragmentTranscript(sequence, coordinate, step, length):
     d = {}
     count = 0
     start = coordinate
     for i in range(0, len(str(sequence)), step):
        v = str(sequence)[(count):(count+length)]
        if len(v) == length:
             k = str(start) + '-' + str(start+len(v)-1)
             d[k] = v
             count += step
             start += step  
        elif len(v) == len(v)%length:
             break 
     return d


def doTranscribeRNA(input_file):
     for record in SeqIO.parse(open(input_file), 'fasta'):
          seq_record = record.seq
          sequence = seq_record.transcribe()
     return sequence

def getEndLocation(input_file):
     records = SeqIO.parse(open(input_file), 'fasta')
     record = next(records)
     parts = record.description.split(':')
     end_location = parts[-2]
     return end_location


def doFastaOutput(input_file, dictionary):
     f = input_file.split('.')
     file_name = f[0] + '_output.' + f[1]
     ofile = open(file_name, 'w')
     for i in dictionary:
          ofile.write('>' + str(i) + '\n' + \
                      str(dictionary[i]) + '\n')
     ofile.close() 
     return 
          

     

