#
# 04/13,10 Reconcile differences between calculon and HPC version.
#


import FileUtility,FastaManager,sys,os,ParseBlast,time,string,math
from time import time

class blast_util:

	##
	# This function is written to automate some aspects of BLAST search and 
	# approval of matching pairs. This involves:
	#  1. read the cluster file generated by ParseBlast.get_qualified.
	#  2. retrieve seqeunces for each cluster, formatdb, blastall
	#  3. parse output and query user interactively, save user input
	#  4. delete temp database, output, and sequence file
	#  5. repeat 2 till no more cluster
	#
	# This is useful in when trying to decide if two things are from the same
	# gene or not semi-manually.
	##	
	def no_name(self,cfile):
		pass
	
	##
	# This is written to very HM clades using outgroup sequences
	#
	# @param clades  tab-delim file with [hs1,...][mm1,....]
	# @param fasta   fasta file with hs and mm sequences
	# @param fasta2  fasta file with the sequences made up the blast db
	# @param all_others are for running blast
	#
	def verify_clade(self,bdir,fasta,fasta2,d,clades,DEBUG):
		
		# read clades into dict: 
		# idx as key, a nested list with [[hs1...][mm1..]] as value
		print "Read clade info..."
		inp = open(clades,"r")
		inl = inp.readline()
		idx = 0
		cdict = {}
		while inl != "":
			inl = self.onl(inl)
			L   = inl.split("\t")
			cdict[idx] = [L[0].split(","),L[1].split(",")]	
			idx += 1	
			inl = inp.readline()
		
		print "Read fasta file..."
		fasta = fm.fasta_to_dict(fasta)
		fasta2= fm.fasta_to_dict(fasta2)
		
		print "Iterate through %i clades:" % len(cdict.keys())
		oup_verify = open(clades+".verify","w")
		for i in cdict:
			print "",i
			# construct TMP.fa and run blast
			oup = open("TMP.fa","w")
			# output hs sequences
			for j in cdict[i][0]:
				oup.write(">%s\n%s\n" % (j,fasta[j][1]))
			# output mm sequences
			for j in cdict[i][1]:
				oup.write(">%s\n%s\n" % (j,fasta[j][1]))
			oup.close()
			# run blast
			self.blast("blastp","TMP.fa",d,"TMP.out",8,bdir,1,F,0,1)
			
			# get top match outgroup sequences
			inp2 = open("TMP.out","r")
			inl2 = inp2.readline()
			subj = {}
			while inl2 != "":
				L = inl2.split("\t")
				if not subj.has_key(L[1]):
					subj[L[1]] = 0
				inl2 = inp2.readline()
			oup = open("TMP2.fa","w")
			for j in subj:
				oup.write(">%s\n%s\n" % (j,fasta2[j][1]))
			oup.close()
			# construct db and run blast again
			os.system("cat TMP.fa TMP2.fa > TMP3.fa")
			os.system("%s/formatdb -i TMP3.fa -n TMP3" % bdir)
			self.blast("blastp","TMP3.fa","TMP3","TMP3.out",8,bdir,1,"F",0,500)
			
			# parse output, get the E VALUEs for 2 things:
			# 1. id of lowest scoring HM pair              -> Xlow  largest  E
			# 2. id of highest scoring H-out or M-out pair -> Xout  smallest E
			# If 1 > 2 : this clade is disqualified
			
			# first construct mdict: id1 id2 as key, id as value
			inp2 = open("TMP3.out","r")
			inl2 = inp2.readline()
			mdict = {}
			while inl2 != "":
				L = inl2.split("\t")
				N = "%s %s" % (L[0],L[1])
				if not mdict.has_key("%s %s" % (L[0],L[1])) and \
				   not mdict.has_key("%s %s" % (L[1],L[0])):
					mdict["%s %s" % (L[0],L[1])] = float(L[2])
				inl2 = inp2.readline()
			
			# get Xlow
			Xlow = ["-","-",100] # [gene1,gene2,transformed evalue]
			                     # transformed evalue is not good, some protein
			                     # are really long and evalue is effectively close
			                     # to zero. Use identity instead.
			# iterate through sp1
			for j in cdict[i][0]:
				# iterate through sp2
				for k in cdict[i][1]:
					if mdict.has_key("%s %s" % (j,k)):
						if Xlow[2] >= mdict["%s %s" % (j,k)]:
							Xlow = [j,k,mdict["%s %s" % (j,k)]]
					elif mdict.has_key("%s %s" % (k,j)):
						if Xlow[2] >= mdict["%s %s" % (k,j)]:
							Xlow = [k,j,mdict["%s %s" % (k,j)]]
					else:
						print "X Score not found:",j,k
			# get Xout
			Xout = ["-","-",0]
			for j in (cdict[i][0]+cdict[i][1]):
				for k in subj:
					if mdict.has_key("%s %s" % (j,k)):
						if Xout[2] < mdict["%s %s" % (j,k)]:
							Xout = [j,k,mdict["%s %s" % (j,k)]]
					elif mdict.has_key("%s %s" % (k,j)):
						if Xout[2] < mdict["%s %s" % (k,j)]:
							Xout = [k,j,mdict["%s %s" % (k,j)]]
					else:
						print "Out score not found:",j,k

			cpass = cfail = 0
			oup_verify.write("%s\t%s\t%s-%s\t%f\t%s-%s\t%f\n" % \
							  (string.joinfields(cdict[i][0],","),
							   string.joinfields(cdict[i][1],","),
							   Xlow[0],Xlow[1],Xlow[2],Xout[0],Xout[1],Xout[2]))
			if DEBUG:
				break
			os.system("rm TMP*")
			
		print "Done!"

	#
	# Take a fasta file, divide it, then call qsub2 to submit jobs
	#
	# @param bdir  blast executable dir
	# @param by    # of files to divide the fasta into
	# @param pm    blast parameters, default is 
	#                p=blastp, m=8, e=1,F=F,v=0,b=2000
	# @param stype sequence type, pep (default) or nt
	# @param btype blast type, default is blastall.
	# @param db    default will be the same as fasta
	# @param pdir  directory where qsub.py is.
	#
	def batch_blast2(self,fasta,bdir,by,pm,stype,btype,db,pdir):
		
		if pm == "":
			pm = "-p blastp -d %s -m 8 -e 1 -F F -v 0 -b 2000" % fasta
		if db == "":
			db = fasta
		
		print "Make database..."
		if os.path.isfile("%s.nin" % db):
			print " already there..." 
		elif stype == "pep":
			os.system("%s/formatdb -i %s" % (bdir,db))
		elif stype == "nt":
			os.system("%s/formatdb -i %s -p F" % (bdir,db))
		else:
			print "Unknown seq type %s, QUIT!" % stype 
		
		print "Divide fasta file..."
		fm.divide(fasta,by)
		
		print "Sumbit jobs..."
		for i in range(1,by+1):
			os.system("python %s/qsub.py -f qsub1 -J batch_blast" % pdir+\
					  		" -c '%s/%s %s -i %s_%i -o %i.out -d %s'"  % \
					  			(bdir,btype,pm,fasta,i,i,db))
		

	#
	# Conduct batch blast search for the fasta files in the wdir
	#
	# Assume blastp,self_db,m=8,bdir=~/bin/blast/,e=1,F=F,v=0,b=2000
	#
	# This also call the best_match function. Coding NOT VERY GOOD.
	#
	# @param bdir blast dir
	# @param wdir sequence dir
	#
	def batch_blast(self,bdir,wdir,pdir,sp1,sp2):
		
		prog = "blastp"
		
		os.system("mkdir _blast")
		os.system("mkdir _csbm")
		
		print "Process fasta file:"
		dlist = os.listdir(wdir)
		countF = 1
		for i in dlist:
			if i[-3:] == ".fa":
				print countF,i 
				countF += 1
				if os.path.isfile("_blast/%s.out" % i):
					if os.path.getsize("_blast/%s.out" % i) != 0:
						print " blast output exists..."
					else:
						print " blast output size is 0"
						os.system("rm %s.out" % i)						
				else:
					self.formatdb(bdir,"%s/%s" % (wdir,i),"T")
					os.system("mv %s/%s.p* ./" % (wdir,i))
					self.blast(prog,"%s/%s" % (wdir,i),i,i+".out",8,bdir,
							   "1","F",0,2000)
					if os.path.getsize("%s.out" % i) != 0:
						os.system("mv %s.out _blast" % i)
					else:
						print "  blast output size is 0"
						os.system("rm %s.out" % i)
					## remove blast db
					os.system("rm %s.p*" % i)
				
				## Call best match
				if os.path.isfile("_csbm/%s.out.cluster" % i):
					if os.path.getsize("_csbm/%s.out.cluster" % i) != 0:
						print " best match done"
					else:
						print " best_match out size is 0"
						os.system("rm %s.out" % i)						
				else:					
					print "python %s/best_match.py " % pdir+\
							"_blast/%s.out %s %s"  % (i,sp1,sp2)
					os.system("python %s/best_match.py " % pdir+\
							  "_blast/%s.out %s %s"  % (i,sp1,sp2))
					if os.path.getsize("_blast/%s.out.cluster" % i) != 0:
						os.system("mv _blast/%s.out.* _csbm"   % i)
					else:
						print "  best_match out size is 0"
						os.system("rm _blast/%s.out.cluster" % i)						
					
		print "Done!"			
	
	#
	# Take a fasta file, divide it, then generate a cmd file for qsub script.
	#
	# @param bdir  blast executable dir
	# @param by    # of files to divide the fasta into
	# @param pm    blast parameters, default is 
	#                p=blastp, m=8, e=1,F=F,v=0,b=2000
	# @param stype database sequence type, pep (default) or nt
	# @param btype blast type, default is blastall.
	# @param db    empty string means fasta should be used.  
	# @param fdir  directory where fasta file is from
	# @param splitfa 1 (yes, default), 0 (no)
	# @param createdb 1 (yes, default), 0 (no)
	#
	def batch_blast_cmd(self,fasta,bdir,by,pm,stype,btype,db,fdir,splitfa,\
				creatdb,outname):
		
		# NOTE!
		print "NOTE: assume qsub location"
		
		if pm == "":
			pm = "-p blastp -d %s -m 8 -e 1 -F F -v 0 -b 2000" % fasta

		if fdir != "":
			if fdir[-1] != "/":
				fdir += "/"
		# Get current working dir if not specified
		else:
			fdir = os.getcwd()+"/"
			print "NOTE: assume working dir:",fdir
				
		if db == "":
			db = fasta
		
		if outname == "":
			outname = fasta
		
		print "Make database..."
		if createdb == 0:
			print " createdb=0, no db created..."
		elif stype == "pep":
			os.system("%s/formatdb -i %s" % (bdir,db))
		elif stype == "nt":
			os.system("%s/formatdb -i %s -p F" % (bdir,db))
		else:
			print "Unknown seq type %s, QUIT!" % stype 
		
		print "Divide fasta file..."
		if splitfa == 0:
			print " splitfa=0, won't divide fasta"
		else:
			fm.divide(fasta,by,0,1)
		
		
		print "Output command lines..."
		oup = open("cmd_%s" % outname, "w")
		for i in range(1,by+1):
			oup.write("%s/%s %s -i %s%s_%i -o %sout.%s_%i -d %s%s\n" % \
					  			(bdir,btype,pm,fdir,fasta,i,fdir,outname,i,fdir,db))
		


	def formatdb(self,d,fasta,p):
		print "   formatdb..."
		os.system("%s/formatdb -i %s -p %s" % (d,fasta,p))
	
	##
	# @param p         blast program
	# @param fasta     input sequence file
	# @param d         database
	# @param o         output file name
	# @param m         output format, set to m=9 because m=8 trucates names.
	# @param bdir blast directory
	# @param e         evalue cutoff, default 1 
	# @param F         filter  
	##	
	def blast(self,p,fasta,d,o,m,bdir,e,F,v,b):
		print "   blast..."
		os.system("%s/blastall -p %s -i %s -d %s -m %i -F %s -e %s -o %s -v %i -b %i" % \
			   	  (bdir,p,fasta,d,m,F,e,o,v,b))

	
	##
	# Batch operation for blast_window. The parameters are mostly the same as 
	# those in blast_window, excepto p, d, o, m which are set to constants.
	#
	# Although no fasta file is asked for, the script NEED fasta files for chr
	# in the working dir.
	#
	# @param block  [block_id][1_geneL][1_geneR][2_geneL][2_geneR]
	##
	def batch_bw(self,block,dat,bdir,w,s,e,n,F):
		
		p     = "tblastx"
		m     = 9
		FLANK = 2500     # add more sequneces to the gene boundaries
		
		# read sequence list into a dict, block id as key, a nested list with:
		#    block_i_a gene name in 0
		#    block_i_b gene name in 1
		bdict = fu.file_to_dict(block,3)
		
		# get the coords
		print "Read block coords..."
		gdict = {}
		for i in bdict:
			for j in bdict[i]:
				if not gdict.has_key(j):
					gdict[j] = []
				
		# dat file: CHR	CODE TYPE ORI L R
		# generate a linked list like dict with chr as key, a dict as value
		# with geneA_L as key, a list as value with 
		# [geneA_R,ORI,NAME,downsream_geneL,downstream_geneR]
		inp = open(dat,"r")
		inl = inp.readline()
		if inl[:3].lower() == "chr":
			inl = inp.readline()
		cdict = {}
		while inl != "":
			L   = inl.split("\t")
			if L[3] == "W":
				L[3] = 1
			else:
				L[3] = -1		
			inl2= inp.readline()
			
			if inl2 != "":
				L2  = inl2.split("\t")
			else:
				break
				
			if gdict.has_key(L[1]):
				gdict[L[1]] = [L[0],int(L[-2]),int(L[-1])]
			if cdict.has_key(L[0]):
				# already has this L, skip
				if cdict[L[0]].has_key(int(L[-2])):
					print " same coord for two entries -",L[-2]
				else:
					# the next gene has the same L, get the L and R for the next
					if int(L[-2]) == int(L2[-2]):
						inl2 = inp.readline()
						if inl2 != "":
							L2 = inl2.split("\t")
						else:
							break
					cdict[L[0]][int(L[-2])] = [int(L[-1]),L[3],L[1],
											   int(L2[-2]),int(L2[-1])]
			else:
				cdict[L[0]] = {int(L[-2]):[int(L[-1]),L[3],L[1],
										   int(L2[-2]),int(L2[-1])]}
			inl = inl2
		
		# iterate through blocks
		print "Iterating blocks.."
		for i in bdict:
			print "",i
			
			if os.path.isfile("%s.out" % i):
				print "  blast done..."
				continue
				
			b1L   = gdict[bdict[i][0]][1]-FLANK
			b1R   = gdict[bdict[i][1]][2]+FLANK
			b2L   = gdict[bdict[i][2]][1]-FLANK
			b2R   = gdict[bdict[i][3]][2]+FLANK
			
			# generate fasta file for block_i_a and b
			# The fasta files has the name chr_id.fa. All should be in the
			# working directory.
			fm.get_stretch2("%s.fa" % i[0],"%i,%i" % (b1L,b1R),0)
			os.system("mv %s.fa.seg.fa TMP1.fa" % i[0])
			fm.get_stretch2("%s.fa" % i[1],"%i,%i" % (b2L,b2R),0)
			os.system("mv %s.fa.seg.fa TMP2.fa" % i[1])
			
			# generate database
			print "Generate database..."
			os.system("%s/formatdb -i TMP2.fa -n TMP -p F" % bdir)
			
			# run blast
			self.blast_window(p,"TMP1.fa","TMP",i+".out",m,bdir,
							   w,s,e,n,F)
			os.system("rm TMP*")
						
			# generate coordinate file for generating graphics
			#   organism \t mm
			#   chr      \t 3
			#   start    \t 1345524323
			#   end      \t 1356664333
			#   gene     \t L \t R \t ORI \t NAME (ori is 1 or -1)
			#   gene     \t ...
			oup_1 = open(i+".coord1","w")
			oup_2 = open(i+".coord2","w")
			oup_1.write("organism\t%s1\nchr\t%s\nstart\t%i\nend\t%i\n" % \
						(i,i[0],b1L,b1R))
			oup_2.write("organism\t%s2\nchr\t%s\nstart\t%i\nend\t%i\n" % \
						(i,i[1],b2L,b2R))
						
			# cdict, key is chr, value is a dict with geneA_L as key, value is
			# [geneA_R,ORI,NAME,downsream_geneL,downstream_geneR]
			# first do block1 genes
			b1L  = b1L + FLANK            # block1 start
			b1R  = b1R - FLANK            # block1 end
			g1L  = gdict[bdict[i][0]][1]  # block1 gene L
			g1R  = cdict[i[0]][g1L][0]    # block1 gene R
			ori1 = cdict[i[0]][g1L][1]    # block1 gene ori
			g1   = cdict[i[0]][g1L][2]    # block1 gene name
			bflag = 1			
			while bflag:
				oup_1.write("gene\t%i\t%i\t%i\t%s\n"%(g1L,g1R,ori1,g1))
				if g1R == b1R:
					bflag = 0	
				g1L  = cdict[i[0]][g1L][3]
				g1R  = cdict[i[0]][g1L][0]
				ori1 = cdict[i[0]][g1L][1]
				g1   = cdict[i[0]][g1L][2]
				
			# do block2 genes
			b2L  = b2L + FLANK            # block2 start
			b2R  = b2R - FLANK            # block2 end
			g2L  = gdict[bdict[i][2]][1]  # block2
			g2R  = cdict[i[1]][g2L][0]
			ori2 = cdict[i[1]][g2L][1]
			g2   = cdict[i[1]][g2L][2]
			bflag = 1			
			while bflag:
				oup_2.write("gene\t%i\t%i\t%i\t%s\n"%(g2L,g2R,ori2,g2))
				if g2R == b2R:
					bflag = 0	
				g2L  = cdict[i[1]][g2L][3]
				g2R  = cdict[i[1]][g2L][0]
				ori2 = cdict[i[1]][g2L][1]
				g2   = cdict[i[1]][g2L][2]
			
			oup_1.close()
			oup_2.close()
			
		print "Done!!"
		
	##
	# @param p         blast program
	# @param i         input sequence file
	# @param d         database
	# @param o         output file name
	# @param m         output format, set to m=9 because m=8 trucates names.
	# @param bdir blast directory    
	# @param w         window size
	# @param s         step size
	# @param e         evalue cutoff, default 1
	# @param n         mega blast, default 1         
	##	
	def blast_window(self,p,fasta,d,o,m,bdir,w,s,e,n,F):
		
		# check if output file is present
		if o == "":
			o = "%s_%i_%s.out" % (fasta,w,s) 
			c = 0
			def check(name,o,c):
				try:
					inp = open(name,"r")
					c += 1
					return [1,o,c]
				except IOError:
					return [0,name,c]
			
			check_again = 1
			while check_again:
				#print check_again,o
				[check_again,o,c] = check("%s.%i" % (o,c),o,c)
		print " output name:",o
		
		# set megablast
		if n == 1:
			n = "T"
		else:
		 	n = "F"
		# set filter
		if F == 1:
			F = "T"
		else:
			F = "F"
		
	
		print "Read seq into dict..."
		seq_dict = fm.fasta_to_dict(fasta,0,1)
		print seq_dict.keys()
		idx      = seq_dict.keys()[0]
		seq      = seq_dict.values()[0]
		print " seq    =",idx
		print " length =",len(seq)
				
		print "Start blast..."
		c = 0
		while c < len(seq):
			#print " %i-%i" % (c+1,c+w)
			segment = seq[c:c+w]
			sname   = "%s_%i-%i" % (idx,c+1,c+w+1)
			oup     = open(sname+".fa","w")
			oup.write(">%s\n%s\n" % (sname,segment))
			
			#print sname
			oup.close()
			# call blast
			#print "%s/blastall -p %s -i %s -d %s -m %i -F F -n T -e %s >> %s" % \
			#	   	  (bdir,p,sname+".fa",d,m,e,o)
			os.system("%s/blastall -p %s -i %s -d %s -m %i -F %s -n T -e %s >> %s" % \
				   	  (bdir,p,sname+".fa",d,m,F,e,o))
			c += s
			
			os.system("rm %s" % (sname+".fa"))
			
		print "Done!"

	
	def bl2seq(self,program,fasta,s1,s2,stdout=1):
	
		fasta = fm.fasta_to_dict(fasta)
		missed = 0
		if not fasta.has_key(s1):
			print "Not found:", s1
			missed = 1
		if not fasta.has_key(s2):
			print "Not found:", s2
			missed = 1
		if missed:
			print "Quit!"
			sys.exit(0)
			
		# generate temp fasta files
		oup = open("TMP.FA1","w")
		oup.write(">%s\n%s\n" % (s1,fasta[s1][1]))
		oup.close()
		oup = open("TMP.FA2","w")
		oup.write(">%s\n%s\n" % (s2,fasta[s2][1]))
		oup.close()
		
		# call bl2seq and parse the output
		if stdout:
			os.system("%s/bl2seq -p %s -i TMP.FA1 -j TMP.FA2 -F F" % \
					   (bdir,program))
		
		else:
			os.system("%s/bl2seq -p %s -i TMP.FA1 -j TMP.FA2 -F F -o %s" % \
					   (bdir,program,"%s_vs_%s.out" % (s1,s2)))   
		os.system("rm TMP*")
	
	#
	# Batch blast sequence pairs and parse output on the fly
	# NOTE: the codes are not very good for running multiple processes at one
	#       time. The tmp files generated will be in the same dir and screw
	#       everything up. So need to deal with this at some point.
	# I don't really understand the implementation before. Don't know why it
	# needs two FASTA files as inputs.
	#
	# I do know. Because this is originally designed for getting the cds coord
	# by cDNA-peptide blast. So input was one cds file, one peptide file.
	# Allow the use of only 1 fasta with all sequences or two as before.
	#
	#      
	def batch_bl2(self, program, pairs, fasta1, fasta2, bdir, outflag=0,
				  wdir="", W=3, top=1, ev=10.0):
		
		t = time()
		if wdir != "" and wdir[-1] != "/":
			wdir += "/"
		elif wdir == "":
			wdir = "tmp_%f/" % t
			os.system("mkdir %s" % wdir)
		
		if ev == "1":
			ev = float(ev)
		
		print "Program  :",program
		print "Pair list:",pairs
		print "Fasta1   :",fasta1
		print "Fasta2   :",fasta2
		print "Blast dir:",bdir
		print "Keep outp:",outflag
		print "Working d:",wdir
		
		# read gene pairs
		print "Read gene pairs..."
		gpairs = fu.file_to_list(pairs,1,"\t")
		print " %i pairs" % len(gpairs)
		
		# see if the list need to be inverted or not
		tlist = []
		for i in gpairs:
			tlist.append([i[0],i[1]])
		gpairs = tlist
			
		# read fasta into dict
		print "Read fasta files..."
		fd1 = fm.fasta_to_dict(fasta1)
		
		has2 = 0
		if fasta2 != "":
			fd2 = fm.fasta_to_dict(fasta2)
			has2 = 1
	
		# see if output need to be saved, if so, generate an empty file
		if outflag:
			oup2 = open("%s.bl2seq.raw" % pairs,"w")
			oup2.close()
		
		# reiterate the list and conduct bl2seq
		print "Do bl2seq..."
		oup = open("%s.bl2seq.out" % pairs,"w")
		oupL= open("%s.log"        % pairs,"w")
		oupL.write("id1\tid2\tmaxL_%ID\tmaxL_alnL\tmaxP_%ID\tmaxP_alnL\n")
		c = 0
		for i in gpairs:
			if c%1e2 == 0:
				print " %i x 100" % (c/100)
			c += 1
			#print i
			
			missed = 0
			if has2:
				try:
					seq1 = fd1[i[0]]
					seq2 = fd2[i[1]]
				except KeyError:
					if not fd1.has_key(i[0]):
						print " %s not in fasta file" % i[0]
					if not fd2.has_key(i[1]):
						print " %s not in fasta file" % i[1]
					missed = 1					
			else:
				try:
					seq1 = fd1[i[0]]
					seq2 = fd1[i[1]]
				except KeyError:
					if not fd1.has_key(i[0]):
						print " %s not in fasta file" % i[0]
					if not fd1.has_key(i[1]):
						print " %s not in fasta file" % i[1]
					missed = 1		
			
			if missed == 1:
				oup.write("%s\t%s\n" % (i[0],i[1]))
				continue
		
			# generate temp fasta files
			oup1 = open("%sTMP.FA1" % wdir,"w")
			oup1.write(">%s\n%s\n" % (i[0],seq1))
			oup1.close()
			oup1 = open("%sTMP.FA2" % wdir,"w")
			oup1.write(">%s\n%s\n" % (i[1],seq2))
			oup1.close()
			
			# call bl2seq and parse the output
			os.system("%s/bl2seq -p %s -i %sTMP.FA1 -j %sTMP.FA2 -F F -o %sTMP.OUT -W %i -e %f"%\
					   (bdir,program,wdir,wdir,wdir,W,ev))
			#print "%s/bl2seq -p %s -i %sTMP.FA1 -j %sTMP.FA2 -F F -e 100 -o %sTMP.OUT"%\
			#		   (bdir,program,wdir,wdir,wdir)
			
			try:
				pb.parse_align2("%sTMP.OUT" % wdir,0)
				olist = fu.file_to_list("%sTMP.OUT.mod" % wdir)
				
				# deal with multiple matches between pairs
				if top == 1:
					topL = 0      # top alignment length
					topP = 0      # top percentage identity
					lineL= ""
					lineP= ""
					for j in olist:
						jlist = j.split("\t")
						# not only it has to be longest, but also a good match
						if int(jlist[3]) > topL:
							topL = int(jlist[3])
							lineL= j
						if int(jlist[2]) > topP:
							topP = int(jlist[2])
							lineP= j
	
					L = lineL.split("\t")
					P = lineP.split("\t")
					if lineL != lineP:
						oupL.write("%s\t%s\t%s\t%s\t%s\t%s\n" % \
											(L[0],L[1],L[2],L[3],P[2],P[3]))
					# write longest one by default
					oup.write("%s\n" % lineL)
				else:
					for j in olist:
						oup.write("%s\n" % j)
				
				if outflag:
					os.system("cat %sTMP.OUT >> %s.bl2seq.raw " % (wdir,pairs))
				
			except TypeError:
				# if blast output do not have any alignment, the pb.parse_algn2
				# will not work because it try to parse the alignment. This
				# happens if the sequence length is too short.
				print " no_blast_out:",i
							
			os.system("rm %sTMP*" % wdir)
		
		os.system("rm -rf tmp_%f/" % t)
		print "Done!"
	
	#
	# This is written for Smith-Waterman alignments between a protein seq and
	# a nucleotide sequence. Mainly for the pseudogene project.
	#
	def batch_sw(self,program,pairs,fasta1,fasta2,bdir,outflag=0,wdir="",ev=""):
		
		if wdir == "":
			wdir = "tmp_%f" % time()
			os.system("mkdir %s" % wdir)
		if wdir[-1] != "/":
			wdir += "/"
		
		#program = "tfasty34"
		flags   = "-A -m 3 -q"
		ev      = float(ev)

		print "Program  :",program
		print "Pair list:",pairs
		print "Fasta1   :",fasta1
		print "Fasta2   :",fasta2
		print "Fasta dir:",bdir
		print "Working d:",wdir
		print "Flags    :",flags
		print "E thres  :",ev

		# read gene pairs
		print "Read gene pairs..."
		gpairs = fu.file_to_list(pairs,1,"\t")
		print " %i pairs" % len(gpairs)
			
		# read fasta into dict
		print "Read fasta files..."
		fd1 = fm.fasta_to_dict(fasta1)
		if fasta2 != "":
			fd2 = fm.fasta_to_dict(fasta2)
			for i in fd2:
				if i not in fd1:
					fd1[i] = fd2[i]
		
		# see if output need to be saved, if so, generate an empty file
		if outflag:
			oup2 = open("%s.sw.raw" % pairs,"w")
			oup2.close()
		
		# reiterate the list and conduct sw alignment
		print "Do sw..."
		oup = open("%s.sw.out" % pairs,"w")
		oup3= open("%s.sw.err" % pairs,"w")
		c = 0
		
		for i in gpairs:
			if len(i) > 2:
				i = i[:2]
			
			if c%1e2 == 0:
				print " %i x 100" % (c/100)
			c += 1
			
			missed = 0
			if not fd1.has_key(i[0]) or not fd1.has_key(i[1]):			
				print " %s not in fasta file" % i
				oup.write("%s\t%s\n" % (i[0],i[1]))
				continue
		
			# generate temp fasta files
			oup1 = open("%sTMP.FA1" % wdir,"w")
			oup1.write(">%s\n%s\n" % (i[0],fd1[i[0]]))
			oup1.close()
			oup1 = open("%sTMP.FA2" % wdir,"w")
			oup1.write(">%s\n%s\n" % (i[1],fd1[i[1]]))
			oup1.close()
			
			# call sw and parse the output
			os.system("%s/%s %s %sTMP.FA1 %sTMP.FA2 > %sTMP.OUT"%\
					   (bdir,program,flags,wdir,wdir,wdir))
			
			# Nested list with:
			#   0      1    2     3   4         5    6      7        8
			#   initn, opt, bits, ev, sw_score, ident, sim, overlap, aln
			results = self.parse_tfasty("%sTMP.OUT" % wdir,1)
			if outflag:
				os.system("cat %sTMP.OUT >> %s.sw.raw " % (wdir,pairs))
			
			# check if there is smith-waterman score, if not, there is no
			# alignment and this sequence will be skipped.
			if results[4] == []:
				print " no alignment:",i
			
			else:
				oup.write("#%s %s\n" % (i[0],i[1]))
				#print i[0],i[1]
				#print results
				# go trough the number of matching fragments
				err = 0
				for j in range(len(results[0])):
					tmp = []
					# go through the returned info
					for k in range(len(results)):
						try:
							tmp.append(results[k][j])
						except IndexError:
							print "problem:",i
							oup3.write("%s\t%s\n" % (i[0],i[1]))
							err = 1
							break
					
					# decide if output should be generated
					if err == 1:
						err = 0
						break
					if float(tmp[3]) <= ev:
						oup.write(">%s\n" % " ".join(tmp[:8]))
						for k in tmp[8]:
							oup.write("%s\n" % k)	
							
			os.system("rm %sTMP*" % wdir)
		
		os.system("rmdir %s" % wdir)
		print "Done!"
		
	#
	# Parse Fasta tfasty format 3 output
	#	
	def parse_tfasty(self,f,call=0):
		inp = open(f)
		inl = inp.readline()
		tag1 = "The best scores are:"
		tag2 = "Smith-Waterman score"
		tag3 = ">>>"
		
		# Nested list with:
		#   0      1    2     3   4         5    6      7        8
		#   initn, opt, bits, ev, sw_score, ident, sim, overlap, aln
		info = [[],[],[],[],[],[],[],[],[]]
				
		foundTag1 = foundTag2 = 0
		seq1L = ""
		seq   = 0
		while inl != "":
			if inl.find(tag3) != -1:
				#print "tag3:",[inl]
				#   1>>>xxxx xxxx aa - xxxx aa
				#            ^^^^
				seq1L = inl[inl.find(">>>"):].split(" ")[1]
			# tag1 in this line
			elif inl.find(tag1) != -1:
				#print "tag1:",[inl]
				foundTag1 = 1
			# tag1 found already
			elif foundTag1:
				# begin of start
				if inl[:-1] != "":
					L = inl[:-1].split(" ")
					# 0         1                 -5  -4    -3  -2   -1
					# seq2_name length(malformed) ori initn opt bits E
					#                  ^^^^^^^^^
					tmp = []
					for j in L:
						if j != "":
							tmp.append(j)
					#print "tmp1:",tmp
					# only take initn and on
					for j in range(-4,0):
						info[j+4].append(tmp[j])
				# the end of score part
				else:
					#print "tag1 end"
					foundTag1 = 0
			elif inl.find(tag2) != -1:
				#print "tag2:",[inl]
				L = inl.split(" ")
				# 0  1  2         3   4  5         6
				# xx xx sw_score; id% xx (similar% (overlap)
				tmp = []
				for j in L:
					if j != "":
						tmp.append(j)
				#print "tmp2:",tmp
				info[4].append(tmp[2][:-1])
				info[5].append(tmp[3][:-1])
				info[6].append(tmp[5][1:-1])
				info[7].append(tmp[-1][1:-2])
				
				# create empty list for alignment
				info[8].append([])	
				# increment foundTag2 for list element counting purpose
				foundTag2 += 1
				
			elif foundTag2 > 0:
				if inl[:-1] != "":
					if inl[0] == ">":
						#print "seq start:",foundTag2
						seq += 1
					elif seq == 1:
						#print "seq1:",inl[:-1]
						if info[8][foundTag2-1] == []:
							info[8][foundTag2-1].append(inl[:-1])
						else:
							info[8][foundTag2-1][0] = \
										info[8][foundTag2-1][0]+inl[:-1]
					elif seq == 2:
						#print "seq2:",inl[:-1]
						if len(info[8][foundTag2-1]) == 1:
							#print ">>"
							info[8][foundTag2-1].append(inl[:-1])
						else:
							#print ">>>"
							info[8][foundTag2-1][1] = \
										info[8][foundTag2-1][1]+inl[:-1]
				else:
					#print "tag2 end"
					# reset sequence count, 0: query, 1: subject
					seq = 0
					
			inl = inp.readline()
		
		if call:
			return info
		else:
			print "\nParsed:"
			for i in info:
				print i
		
	
	#
	# Adjust the sequence coordinates based on gap coordinates
	#
	# @param gap   the ".gap" file generated by ParseBlast.pares_gap
	# @param pair  sequence id pairs. The id should followed by "|dom_name"
	# @param coord [seqid][L][R][dom]. The sequence id here is plain.
	#	
	def verify_match(self,gap,pair,coord):
		
		# read domain coordinates
		inp = open(coord,"r")
		inl = inp.readline()
		ddict = {} # dom as key, value a dict w/ gene_id as key, value [L,R]
		while inl != "":
			L = self.onl(inl).split("\t")
			if ddict.has_key(L[-1]):
				# multiple domain thing, NOT taken care of.
				if ddict[L[-1]].has_key(L[0]):
					ddict[L[-1]][L[0]] = "-"
				else:
					ddict[L[-1]][L[0]] = [int(L[1]),int(L[2])]
			else:
				ddict[L[-1]] = {L[0]:[int(L[1]),int(L[2])]}	
			inl = inp.readline()
		
		inp = open(pair,"r")
		
	
	
	
	def onl(self,line):
		if line[-2:] == "\r\n":
			line = line[:-2]
		elif line[-1] == "\n":
			line = line[:-1]
		return line

	
	def help(self):
		print " -f function"	
		print "     batch_bw - batch operation of blast_window. REQUIRES: b,"
		print "        dat. OPTIONAL: bdir,w,s,e,n"
		print "     batch_blast - require dir with sequence file, REQUIRES: D,"
		print "        s1, s2, OPTIONAL: bdir,pdir"
		print "     batch_blast2 - divide a fasta and run multiple blast. NEED:"
		print "        fasta, OPT: bdir, by, pm, stype, btype,db, pdir"
		print "     batch_bl2 - batch blast a list of pairs. REQUIRES: p,g,i"
		print "        OPTIONAL: j,bdir,d,wdir,W,top,e"
		print "     batch_blast_cmd - divide a fasta and generate a cmd list." 
		print "        NEED:fasta, OPT: bdir, by, pm, stype, btype,db, fdir"
		print "        splitfa, createdb, outname"
		print "     batch_sw - do batch tfasty with smith-waterman algorithm"
		print "        NEED: p,g,i, OPT: j, bdir, d, e"
		print "     bl2seq - take two sequences from a fasta file and do bl2seq"
		print "        REQUIRES: p,fasta,s1,s2"
		print "     blast_window - blast windows of a given sequence. REQUIRES:"
		print "        p,i,d. OPTIONAL: o,m,bdir,w,s,e,n"
		print "     parse_tfasty - NEED: infile"
		print "     verify_clade - a special function for blasting each"
		print "        ancestral units defined against database made of seq"
		print "        from outgroup organsisms and verify if within clade"
		print "        matches are better than to outgroups. REQUIRES:bdir"
		print "        i,j,c,d"
		print "     verify_match - NEED: gap,g,c"
		print " -p  blast program name"
		print " -createdb 1 (yes, default), 0 (no)"
		print " -g  pair list. For verify_match, check the code doc."
		print " -i  first fasta file, should be the index 0 of the pairs." 
		print "     For batch_bl2 this file can have all sequences if the run" 
		print "        is not about getting protein sequences out of nt"
		print " -j  second fasta file, should be the index 1 of the pairs"
		print " -r  which of each pair to use as query, 0 (default) or 1" 
		print " -fasta  query seq file"
		print " -d  subject database, for batch_bl2, this is the flag for"
		print "     keeping the bl2seq or tfasty output (1) or not (0, default)"
		print " -m  ouput format, default 8"
		print " -bdir  blast folder directory, default /usr/bin/blast"
		print " -pdir  python code dir, default ~/codes/"
		print " -fdir  fasta file folder where database and outputs will be."
		print " -w  window size, default 1000"
		print " -W  word size, default 3"
		print " -by divide fasta file by, default 1"
		print " -pm parameter, default: '-p blastp -m 8 -e 1 -F F -v 0 -b 2000"
		print " -outname output postfix for batch_blast_cmd"
		print " -s  step size, default 1000"
		print " -splitfa 1 (divide fasta, default), 0 (no)"
		print " -stype database sequence type, pep (default) or nt"
		print " -btype blast type, default blastall."
		print " -db database or fasta for database. If not passed, assume it is"
		print "     a self blast search"
		print " -e  evalue cutoff, default 1"
		print " -n  mega blast, default 1 (which is true)"
		print " -F  low-complexity filter, default 0 (false)"
		print " -s1 sequence name 1, or species abbr 1 for batch_blast"
		print " -s2 sequence name 2, or species abbr 2"
		print " -b  block infor with [gene][ori][L][R][block_id]"
		print " -dat dat file generated by DatParser.parse_coord"
		print " -c  clade file with [o1g1,o1g2....][o2g1,o2g2...]. For"
		print "        verify_match, this is the coord file"
		print " -D  directory with sequence files"
		print " -gap Gap file generated by ParseBlast.parse_gap"
		print " -wdir working dir, default ./"
		print " -infile tfasty output"
		print " -top return top match only, default 1 (true)"
		print ""
		sys.exit(0)


if __name__ == '__main__':
	
	util = blast_util()
	fm   = FastaManager.fasta_manager()
	fu   = FileUtility.file_util()
	pb   = ParseBlast.parser()
	f    = fasta = p = o = g = fasta2 = s1 = s2 = b = dat = c = gap = \
		   D = wdir = infile = pm = db = fdir = outname = ""
	d    = 0
	m    = 9
	w    = 5000
	W    = 3
	s    = w
	e	 = "1"
	n	 = 1
	r    = 0
	F    = 0
	DEBUG = 0
	by   = 1
	bdir = "/usr/bin/blast"
	pdir = "~/codes/"
	stype = "pep"
	btype = "blastall"
	top  = 1
	splitfa = 1
	createdb = 1
	
	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-f":
			f   = sys.argv[i+1]
		elif sys.argv[i] == "-p":
			p   = sys.argv[i+1]
		elif sys.argv[i] == "-fasta":
			fasta   = sys.argv[i+1]
		elif sys.argv[i] == "-d":
			d   = int(sys.argv[i+1])
		elif sys.argv[i] == "-o":
			o   = sys.argv[i+1]
		elif sys.argv[i] == "-m":
			m   = int(sys.argv[i+1])
		elif sys.argv[i] == "-w":
			w   = int(sys.argv[i+1])
		elif sys.argv[i] == "-s":
			s   = int(sys.argv[i+1])
		elif sys.argv[i] == "-bdir":
			bdir = sys.argv[i+1]
		elif sys.argv[i] == "-pdir":
			pdir = sys.argv[i+1]
		elif sys.argv[i] == "-fdir":
			fdir       = sys.argv[i+1]
		elif sys.argv[i] == "-e":
			e   = float(sys.argv[i+1])
		elif sys.argv[i] == "-n":
			n   = int(sys.argv[i+1])
		elif sys.argv[i] == "-i":
			fasta    = sys.argv[i+1]
		elif sys.argv[i] == "-j":
			fasta2    = sys.argv[i+1]
		elif sys.argv[i] == "-g":
			g         = sys.argv[i+1]
		elif sys.argv[i] == "-r":
			r         = sys.argv[i+1]
		elif sys.argv[i] == "-s1":
			s1        = sys.argv[i+1]
		elif sys.argv[i] == "-s2":
			s2        = sys.argv[i+1]
		elif sys.argv[i] == "-b":
			b        = sys.argv[i+1]
		elif sys.argv[i] == "-dat":
			dat        = sys.argv[i+1]
		elif sys.argv[i] == "-c":
			c          = sys.argv[i+1]
		elif sys.argv[i] == "-DEBUG":
			DEBUG      = sys.argv[i+1]
		elif sys.argv[i] == "-D":
			D          = sys.argv[i+1]
		elif sys.argv[i] == "-gap":
			gap        = sys.argv[i+1]
		elif sys.argv[i] == "-wdir":
			wdir       = sys.argv[i+1]
		elif sys.argv[i] == "-infile":
			infile     = sys.argv[i+1]
		elif sys.argv[i] == "-by":
			by         = int(sys.argv[i+1])
		elif sys.argv[i] == "-pm":
			pm         = sys.argv[i+1]
		elif sys.argv[i] == "-stype":
			stype         = sys.argv[i+1]
		elif sys.argv[i] == "-btype":
			btype         = sys.argv[i+1]
		elif sys.argv[i] == "-db":
			db         = sys.argv[i+1]
		elif sys.argv[i] == "-W":
			W         = int(sys.argv[i+1])
		elif sys.argv[i] == "-top":
			top         = int(sys.argv[i+1])
		elif sys.argv[i] == "-splitfa":
			splitfa     = int(sys.argv[i+1])
		elif sys.argv[i] == "-createdb":
			createdb    = int(sys.argv[i+1])
		elif sys.argv[i] == "-outname":
			outname    = sys.argv[i+1]
		else:
			print "Unknown flag:",sys.argv[i]
			util.help()
			
	if f == "blast_window":
		if p == "" or fasta1 == "" or d == "":
			print "\nNeed blast program name, input sequence file, and database\n"
			util.help()
		util.blast_window(p,fasta1,d,o,m,bdir,w,s,e,n,F)
	elif f == "batch_bw":
		if b == "" or dat == "":
			print "\nNeed block, dat files\n"
			util.help()
		util.batch_bw(b,dat,bdir,w,s,e,n,F)
	elif f == "batch_bl2":
		if "" in [p,g,fasta]:
			print "\nRequires program, gene pairs, fasta"
			util.help()
		util.batch_bl2(p,g,fasta,fasta2,bdir,d,wdir,W,top,e)
	elif f == "batch_sw":
		if "" in [p,g,fasta]:
			print "\nRequires program, gene pairs, fasta"
			util.help()
		util.batch_sw(p,g,fasta,fasta2,bdir,d,wdir,e)
	elif f == "bl2seq":
		if "" in [p,fasta,s1,s2]:
			print "\nRequires program, fasta, ane two sequence names"
			util.help()
		util.bl2seq(p,fasta,s1,s2)
	elif f == "verify_clade":
		if "" in [bdir,fasta1,fasta2,d,c]:
			print "\nRequires bdir, 2 fasta, clades, and database name"
			util.help()
		util.verify_clade(bdir,fasta1,fasta2,d,c,DEBUG)		
	elif f == "batch_blast":
		if "" in [D,s1,s2]:
			print "\nRequires sequence dir and species abbrv."
			util.help()
		util.batch_blast(bdir,D,pdir,s1,s2)	
	elif f == "batch_blast2":
		if "" in [fasta]:
			print "\nRequires fasta file"
			util.help()
		util.batch_blast2(fasta,bdir,by,pm,stype,btype,db,pdir)	
	elif f == "batch_blast_cmd":
		if "" in [fasta]:
			print "\nRequires fasta file"
			util.help()
		util.batch_blast_cmd(fasta,bdir,by,pm,stype,btype,db,fdir,splitfa,
							createdb,outname)	
	elif f == "verify_match":
		if "" in [gap,g,c]:
			print "\nRequires gap, pair, and coord files"
			util.help()
		util.verify_match(gap,g,c)
	elif f == "parse_tfasty":
		if infile == "":
			print "\nNeed tfasty output"
		util.parse_tfasty(infile)
	else:
		print "\nUnknown function...\n"
		util.help()
	
	
	
