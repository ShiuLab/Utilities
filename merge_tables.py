"""
Script to merge data frames using an ID column

python merge_tables.py [HOW] [DataFile1] [Col_ID1] [DataFile2] [Col_ID2] etc...

HOW: 
left: use only keys from left frame, similar to a SQL left outer join; preserve key order
right: use only keys from right frame, similar to a SQL right outer join; preserve key order
outer: use union of keys from both frames, similar to a SQL full outer join; sort keys lexicographically
inner: use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys

If the data frame has a header row, use the name of the column as the Col_ID
If the data frame does not have a header use the column number (starting with 0) as the Col_ID

"""

import os, sys
import pandas as pd

if len(sys.argv) <= 1:
    print(__doc__)
    exit()

dfs = []
on = []

# Read in arguments to populate dfs and (merge) on lists
HOW = sys.argv[1]
SAVE = sys.argv[2]

for i in range (3,len(sys.argv),2):
  dfs.append(sys.argv[i])
  try:
    on.append(int(sys.argv[i+1]))
  except:
    on.append(sys.argv[i+1])


print("Merging %i dataframes using %s join." % (len(dfs), HOW))

count = 1
for i in range(len(dfs)):  
  if isinstance(on[i], int):
    d_temp = pd.read_csv(dfs[i], header=None, sep=None, engine='python')
  else:
    d_temp = pd.read_csv(dfs[i], header=0, sep=None, engine='python')

  
  if count == 1:
    d = d_temp.copy()
  
  else:
    d = d.merge(right=d_temp, how=HOW, left_on=on[0], right_on=on[i])

  count += 1

print(d.head())

d.to_csv(SAVE, header=True, index=False)

