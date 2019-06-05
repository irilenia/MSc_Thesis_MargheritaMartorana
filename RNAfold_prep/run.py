from functions import *

# Set Parameters
fragment_step = 1
fragment_length = 20

file = 'gene.fasta'

# Run
rna_seq = doTranscribeRNA(file)
start_location = getStartLocation(file)
fragments_dictionary = doFragmentTranscript(\
     rna_seq, start_location, fragment_step, fragment_length)
doFastaOutput(file, fragments_dictionary)




