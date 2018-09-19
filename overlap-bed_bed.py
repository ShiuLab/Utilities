
#import modules
import sys
import OverlapFunctions as ov

def print_help():
	print'''
inp1 = reference BED file
inp2 = target BED file
inp3 = (optional) string for output
'''

def main():
	if len(sys.argv) == 1 or "-h" in sys.argv:
		print_help()
		sys.exit()
	
	try:
		ref_bed_file = sys.argv[1]
		tar_bed_file = sys.argv[2]
		if len(sys.argv) == 4:
			output_str = sys.argv[3]
		else:
			output_str = tar_bed_file.split("/")[-1]
	except:
		print_help()
		print "Error reading arguments, quitting!"
		sys.exit()
	
	arg_dict = {"Ref":[ref_bed_file,0,1,2], "Tar":[tar_bed_file,0,1,2], "Span": 10000, "LnOut": "%s.%s_overlap.lines"%(ref_bed_file,output_str), "NumOut": "%s.%s_overlap.num"%(ref_bed_file,output_str)}
	arg_dict["Ref"][0] = [ln.strip() for ln in open(arg_dict["Ref"][0],"r").readlines()] 
	arg_dict["Tar"][0] = [ln.strip() for ln in open(arg_dict["Tar"][0],"r").readlines()]
	[overlap_dict,sorted_keys] = ov.FindOverlapsByIndexing(arg_dict["Ref"],arg_dict["Tar"],arg_dict["Span"])
	ov.WriteOverlapLines(overlap_dict,sorted_keys,arg_dict["LnOut"])
	ov.WriteOverlapNumber(overlap_dict,sorted_keys,arg_dict["NumOut"])

if __name__ == "__main__":
	main()
