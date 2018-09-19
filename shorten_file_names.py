
#import modules
import sys
import os

def print_help():
	print'''
inp1 = directory with files
inp2 = find/replace file
'''

def make_rename_dict(fl):
	d = {}
	inp = open(fl)
	for line in inp:
		if not line.startswith("#"):
			print line
			find,replace = line.split("\t")
			d[find] = replace.strip()
	inp.close()
	return d

def rename_files(dr,rnm_d):
	fp_dr = os.path.abspath(dr)
	for fl in os.listdir(dr):
		fp_fl = fp_dr+"/"+fl
		renamed = False
		new_nm = fp_fl
		for old_rn in rnm_d.keys():
			if old_rn in fl:
				new_rn = rnm_d[old_rn]
				new_nm = new_nm.replace(old_rn,new_rn)
				renamed = True
		if renamed == True:
			print
			print "OLD:",fp_fl.split("/")[-1]
			print "NEW:",new_nm.split("/")[-1]
			os.system("mv %s %s"%(fp_fl,new_nm))

def main():
	if len(sys.argv) == 1 or "-h" in sys.argv:
		print_help()
		sys.exit()
	
	try:
		in_dir = sys.argv[1]
		find_repl_file = sys.argv[2]
	except:
		print_help()
		print "Error reading arguments, quitting!"
		sys.exit()
	
	rename_dict = make_rename_dict(find_repl_file)
	print rename_dict
	rename_files(in_dir,rename_dict)

if __name__ == "__main__":
	main()
