
#import modules
import sys
import fn

def print_help():
	print'''
inp1 = bed-bed overlap file
'''

def calc_percent_overlap(fl):
	d = {}
	inp = open(fl)
	out = open(fl+".per_over","w")
	for line in inp:
		if not line.startswith("#"):
			lineLst = line.strip().split("\t")
			sq = lineLst[0]
			st1 = int(lineLst[1])
			en1 = int(lineLst[2])
			loc_id = "%s|%s-%s"%(sq,st1,en1)
			if lineLst[4] != "NA":
				st2 = int(lineLst[4])
				en2 = int(lineLst[5])
				per_ovr1,cov_start1,cov_end1 = fn.calc_percent_coverage(st1,en1,st2,en2)
				per_ovr2,cov_start2,cov_end2 = fn.calc_percent_coverage(st2,en2,st1,en1)
			else:
				per_ovr1 = "NA"
				per_ovr2 = "NA"
			new_line = "%s\t%s\t%s\n"%(line.strip(),per_ovr1,per_ovr2)
			out.write(new_line)
			# if loc_id not in d:
				# d[loc_id] = [per_ovr1,new_line]
				# if loc_id == "Zmay_Chr2|19134109-19138815":
					# print d[loc_id]
			# else:
				# prev_per_ovr = d[loc_id][0]
				# if prev_per_ovr == "NA":
					# pass
				# else:
					# if per_ovr1 > prev_per_ovr:
						# d[loc_id] = [per_ovr1,new_line]
	out.close()
	inp.close()
	
	# out = open(fl+".per_over","w")
	# for loc_id in d:
		# out.write(d[loc_id][1])
	# out.close()

def main():
	if len(sys.argv) == 1 or "-h" in sys.argv:
		print_help()
		sys.exit()
	
	try:
		overlap_file = sys.argv[1]
	except:
		print_help()
		print "Error reading arguments, quitting!"
		sys.exit()
	
	calc_percent_overlap(overlap_file)

if __name__ == "__main__":
	main()
