
"""
PURPOSE:
Use -drop or -keep to drop or keep the column names in that list from the dataframe (-df)

INPUTS:
  
  REQUIRED:
  -df       Dataframe
  -keep or -drop    File with list of column names to keep or drop

  OPTIONAL:
  -SAVE   Name to save output as. Default = df_keep or df_drop
  -ID     Name of ID column to keep that isn't in the -keep file
  -sep    Seperator (Default = tab)

"""
import pandas as pd
import os, sys

KEEP = 'f'
DROP = 'f'
SAVE = 'f'
SEP = '\t'
ID = 'f'

for i in range (1,len(sys.argv),2):
  if sys.argv[i].lower() == "-df":
    DF = sys.argv[i+1]
  if sys.argv[i].lower() == "-keep":
    KEEP = sys.argv[i+1]
  if sys.argv[i].lower() == "-drop":
    DROP = sys.argv[i+1]
  if sys.argv[i].lower() == "-save":
    SAVE = sys.argv[i+1]
  if sys.argv[i].lower() == "-sep":
    SEP = sys.argv[i+1]
  if sys.argv[i].lower() == "-id":
    ID = sys.argv[i+1]

if len(sys.argv) <= 1:
  print(__doc__)
  exit()

if SAVE != 'f':
  SAVE_NAME = SAVE
else:
  if KEEP != 'f':
    SAVE_NAME = DF + "_" + KEEP
  elif DROP != 'f':
    SAVE_NAME = DF + '_' + DROP
  else:
    print('Need to specifiy column names to keep or drop!')
    exit()

print('Reading in dataframe...')
d = pd.read_csv(DF, header=0, sep=SEP)

if KEEP != 'f':
  with open(KEEP) as f:
    keep = f.read().splitlines()
  if ID != 'f':
    d = d[[ID] + keep]
  else:
    d = d[keep]

if DROP != 'f':
  with open(DROP) as f:
    drop = f.read().splitlines()   

  d = d.drop(DROP, axis=1)

print(d.head())
d.to_csv(SAVE_NAME, sep=SEP, index=False)
print('Finished!')
