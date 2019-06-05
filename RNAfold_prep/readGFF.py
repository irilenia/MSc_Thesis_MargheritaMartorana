import gffutils
import pandas as pd

# GFF file to use, needs long path
fn = gffutils.example_filename('/Users/Marghina/Desktop/MSc/Project/scripts/gene.gff')

# Create the database
db = gffutils.create_db(fn, dbfn='dcx.db', force=True, \
                        keep_order=True, merge_strategy='merge',\
                        sort_attribute_values=True)
