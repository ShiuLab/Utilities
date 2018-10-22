"""
PURPOSE:
Generate Venn Diagrams with number of overlapping strings from different lists


INPUT:
  -files       Comma separated list of files with lists to compare
  -ids         Comma separated list of what to call each file (same order as scores)
  -save        Save name

OUTPUT:
  -save_venn.txt           Venn diagram
  -save_ovlp_matrix.txt    Matrix with counts for overlaps between _scores files


AUTHOR: Christina Azodi

python ~/GitHub/Utilities/plot_venn.py -files featsel_YLD_EN_50,featsel_YLD_RF_50,featsel_YLD_BA_50 -ids EN,RF,BayA -save plot_venn_FeatSelp50

"""
import pandas as pd
import numpy as np
from collections import defaultdict
import sys, os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import venn

save = 'plot_venn'

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-files":
    FILES = sys.argv[i+1]
  if sys.argv[i] == "-ids":
    IDS = sys.argv[i+1]
  if sys.argv[i] == '-save':
    save = sys.argv[i+1]
if len(sys.argv) <= 1:
  print(__doc__)
  exit()

files = FILES.strip().split(',')
n_comparing = len(files)
ids = IDS.strip().split(',')

# Read all files into list of lists with ID as key
comp = defaultdict(list)
count = 0
for filex in files:
  name = ids[count]
  with open(filex) as f:
    file_list = f.read().splitlines()
  comp[name] = file_list
  count += 1


df = pd.DataFrame(0, index=ids, columns=ids)

# Make a table of overlaps
ids2 = ids[:]
for A in ids:
  for B in ids2:
    if A == B:
      df[A][B] = len(comp[A])
    else:
      df[A][B] = len(list(set(comp[A]).intersection(comp[B])))
  ids2.remove(A)
print(df)

#df.to_csv(out1)


### Output a list of overlaps for making a venn diagram in R "VennDiagram"
def intersection(*listas):
  return set(listas[0]).intersection(*listas[1:])

if n_comparing == 2:
  
  n12 = len(intersection(comp[ids[0]], comp[ids[1]]))
  genes2 =  (comp[ids[1]])
  area2 = len(genes2) 
  genes1 =  (comp[ids[0]])
  area1 = len(genes1)
  unique_2 = list(set(genes2) - set(genes1))
  unique_1 = list(set(genes1) - set(genes2))
  out1.write('\n\nGenes only correctly predicted as positive by:\n')
  out1.write('%s: %s\n%s: %s' % (str(ids[0]), ','.join(str(x) for x in unique_1), ids[1],','.join(str(x) for x in unique_2))) 
  venn_diag_list = [area1, area2, n12]

if n_comparing == 3:
  n123 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]]))
  n12 = len(intersection(comp[ids[0]], comp[ids[1]]))
  n13 = len(intersection(comp[ids[0]], comp[ids[2]])) 
  n23 = len(intersection(comp[ids[1]], comp[ids[2]]))  
  area3 = len(comp[ids[2]]) 
  area2 = len(comp[ids[1]]) 
  area1 = len(comp[ids[0]])  
  venn_diag_list = [area1, area2, area3, n12, n13, n23, n123]


if n_comparing == 4:
  n1234 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]]))
  n234 = len(intersection(comp[ids[1]], comp[ids[2]], comp[ids[3]]))
  n134 = len(intersection(comp[ids[0]], comp[ids[2]], comp[ids[3]]))
  n124 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[3]]))
  n123 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]]))
  n12 = len(intersection(comp[ids[0]], comp[ids[1]]))
  n13 = len(intersection(comp[ids[0]], comp[ids[2]])) 
  n14 = len(intersection(comp[ids[0]], comp[ids[3]])) 
  n23 = len(intersection(comp[ids[1]], comp[ids[2]])) 
  n24 = len(intersection(comp[ids[1]], comp[ids[3]])) 
  n34 = len(intersection(comp[ids[2]], comp[ids[3]])) 
  area4 = len(comp[ids[3]])
  area3 = len(comp[ids[2]]) 
  area2 = len(comp[ids[1]]) 
  area1 = len(comp[ids[0]])  
  venn_diag_list = [area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234]

if n_comparing == 5:
  n12345 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]], comp[ids[4]]))
  n2345 = len(intersection(comp[ids[1]], comp[ids[2]], comp[ids[3]], comp[ids[4]]))
  n1345 = len(intersection(comp[ids[0]], comp[ids[2]], comp[ids[3]], comp[ids[4]]))
  n1245 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[3]], comp[ids[4]]))
  n1235 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[4]]))
  n1234 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]]))
  n345 = len(intersection(comp[ids[2]], comp[ids[3]], comp[ids[4]]))
  n245 = len(intersection(comp[ids[1]], comp[ids[3]], comp[ids[4]]))
  n235 = len(intersection(comp[ids[1]], comp[ids[2]], comp[ids[4]]))
  n234 = len(intersection(comp[ids[1]], comp[ids[2]], comp[ids[3]]))
  n145 = len(intersection(comp[ids[0]], comp[ids[3]], comp[ids[4]]))
  n135 = len(intersection(comp[ids[0]], comp[ids[2]], comp[ids[4]]))
  n134 = len(intersection(comp[ids[0]], comp[ids[2]], comp[ids[3]]))
  n125 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[4]]))
  n124 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[3]]))
  n123 = len(intersection(comp[ids[0]], comp[ids[1]], comp[ids[2]]))
  n12 = len(intersection(comp[ids[0]], comp[ids[1]]))
  n13 = len(intersection(comp[ids[0]], comp[ids[2]])) 
  n14 = len(intersection(comp[ids[0]], comp[ids[3]])) 
  n15 = len(intersection(comp[ids[0]], comp[ids[4]])) 
  n23 = len(intersection(comp[ids[1]], comp[ids[2]])) 
  n24 = len(intersection(comp[ids[1]], comp[ids[3]])) 
  n25 = len(intersection(comp[ids[1]], comp[ids[4]])) 
  n34 = len(intersection(comp[ids[2]], comp[ids[3]])) 
  n35 = len(intersection(comp[ids[2]], comp[ids[4]]))
  n45 = len(intersection(comp[ids[3]], comp[ids[4]])) 
  area5 = len(comp[ids[4]]) 
  area4 = len(comp[ids[3]])
  area3 = len(comp[ids[2]]) 
  area2 = len(comp[ids[1]]) 
  area1 = len(comp[ids[0]])  
  venn_diag_list = [area1, area2, area3, area4, area5, n12, n13, n14, n15,n23, n24, n25,
   n34, n35, n45, n123, n124, n125, n134, n135, n145, n234, n235, n245, n345, n1234, 
   n1235, n1245, n1345, n2345, n12345]

print('Overlap list needed for VennDiagram in R')
print(venn_diag_list)
#out1.write('\n\n## Overlap list needed for VennDiagram in R:\n\n%s' % venn_diag_list)

if n_comparing == 5:
  labels = venn.get_labels([comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]], comp[ids[4]]], fill = ['number'])
  fig, ax = venn.venn5(labels, names = ids)
elif n_comparing == 4:
  labels = venn.get_labels([comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]]], fill = ['number'])
  fig, ax = venn.venn4(labels, names = ids)
elif n_comparing == 3:
  labels = venn.get_labels([comp[ids[0]], comp[ids[1]], comp[ids[2]]], fill = ['number'])
  fig, ax = venn.venn3(labels, names = ids)
elif n_comparing == 2:
  labels = venn.get_labels([comp[ids[0]], comp[ids[1]]], fill = ['number'])
  fig, ax = venn.venn2(labels, names = ids)
filename = save+'_pred_compared.pdf'
fig.savefig(filename)