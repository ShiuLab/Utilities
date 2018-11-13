
import sys, os, time
from random import randint

class qsub_hpc:

    def write_sh(self,cmd,jobname,sidx,p,h,m,mem,email,wd,mo,pre,post,a,inta,nnode,ngpu,array,devnode):
        # Write header
        oup = open("%s%i.sb" % (jobname,sidx),"w")
        oup.write('#!/bin/bash')
        oup.write('\n########## Define Resources Needed with SBATCH Lines ##########\n\n')
        oup.write("#SBATCH -n %i -c %i --gres=gpus:%i\n" % (nnode, p, ngpu))
        oup.write("#SBATCH --time=%i:%i\n" % (h,m))
        oup.write("#SBATCH --mem=%iG\n" % (mem))
        if email != "":
            oup.write("#SBATCH --mail-user=%s\n" % email)
            oup.write("#SBATCH --mail-type=FAIL,BEGIN,END\n")
        if jobname != "":
            oup.write("#SBATCH -J %s\n" % jobname)
        if a != "":
            oup.write("#SBATCH -A %s\n" % a)
        if inta != "":
            oup.write("srun --pty /bin/bash")
        if array != "":
            oup.write("#SBATCH --array=%s\n" % array)
        if devnode != "":
            oup.write("#SBATCH -C %s\n" % devnode)
        
        oup.write('\n########## Command Lines to Run ##########\n\n')

        # Set working directory
        if wd != "":
            oup.write("cd %s\n" % wd)

        # Load modules
        if mo != []:
            for j in mo:
                oup.write("module load %s\n" % j)
                
        # Write cmd lines defined in the pre parameter
        if pre != "":
            oup.write("%s\n" % "".join(open(pre).readlines()))
        
        # Write main cmd line
        oup.write("%s\n\n" % cmd.strip())

        # Write cmd lines defined in the post parameter
        if post != "":
            oup.write("%s\n" % "".join(open(post).readlines()))
        oup.close()
    
    def submit(self,jobs,sidx,wtime,mem,jobname,p,email,logdir,a,inta,nnode,ngpu,array,devnode,wdir="",
              module="",pre="",post=""):
        runtype = 0
        if type(jobs) != list:
            runtype = 1
            inp = open(jobs)
            inl = inp.readlines()
            tmp = []
            for i in inl:
                i = i.strip()
                tmp.append(i)
            jobs = tmp

        if logdir == 0:
            logdir = "-o tmp.o -e tmp.e "
        else:
            logdir = ""
    
        # read modules to be loaded
        if module == "":
            mo = []
        elif os.path.isfile(module):
            inp = open(module)
            mo  = inp.readlines()
        else:
            mo  = module.split(',')

        h = wtime/60
        m = wtime%60
    
        for i in jobs:
            if i.strip() != '' and i[0] != "#":
                print("  job %i" % sidx)
                self.write_sh(i,jobname,sidx,p,h,m,mem,email,wdir,mo,pre,post,a,inta,nnode,ngpu,array,devnode)
                os.system("chmod 755 %s%i.sb" % (jobname,sidx))
                os.system("sbatch %s%s%i.sb"    % (logdir,jobname,sidx))  
                sidx += 1
    
    def queue(self,jcommand,stime,nsub,juser,jrange,wtime,mem,jobname,p,email,
            logdir,a,inta,nnode,ngpu,array,devnode,wdir="",module="",pre="",post=""):
        jdict = {}
        inp = open(jcommand)
        jobs = inp.readlines()
        if jrange != "*":
            jrange = jrange.split("-")
            jobs  = jobs[int(jrange[0])-1:int(jrange[1])]
        
        oup = open("%s.log" % jcommand,"w")                 
        j = 0
        while j < len(jobs): 
            oup.write("%s\t" % time.ctime())
            print("%s\t" % time.ctime())
            
            #####
            # Check number of jobs in quene so far
            #####
            rint = randint(1,1e6)
            os.system("squeue -l -u %s > TMP.%i" % (juser,rint))
            time.sleep(2)
            inp = open("TMP.%i" % rint)
            inl = inp.readlines()
            if len(inl) != 0:
                currj = len(inl)-2
            else:
                currj = 0
            inp.close()
            os.system("rm TMP.%i" % rint)   
            #####
            # Submit more jobs so there is always nsub number of jobs in quene
            #####
            if currj < nsub:
                if array != '':
                    array_start, array_stop = array.strip().split('-')
                    array_size = int(array_stop) - int(array_start) + 1
                    jseg = jobs[j:j+int((nsub-currj) / array_size)]
                else:
                    jseg = jobs[j:j+(nsub-currj)]

                oup.write("submit %i\t" % len(jseg))
                print("current %i, submit %i\n" % (currj,len(jseg)))
                self.submit(jseg,j+1,wtime,mem,jobname,p,email,logdir,a,inta,nnode,ngpu,array,devnode,wdir,module)
                j += len(jseg)

            else:
                oup.write("waiting\t")
                print("waiting\t")
            oup.write("submited so far: %i\n" % j)
            print("submited so far: %i\n" % j)
            time.sleep(stime)
        #submit(self,jobs,sidx,wtime,mem,jobname,p,email,logdir,a,i,nnode,ngpu,array,wdir="",module="",pre="",post=""):
    print("Done!")

    def qdel(self,duser,drange):
        #if "" not in [duser,drange]:
        #   print("Specify either user OR range! Quit!...\n"
        #   sys.exit(0)
        
        print("User :",duser)
        print("Range:",drange)
        
        if duser != "":
            os.system("squeue -l -u %s > TMP_qdel" % duser)
            inp = open("TMP_qdel")
            inl = inp.readlines()[5:]
            dlist = []
            for i in inl:
                dlist.append(i.split(".")[0])
            print("Kill:")
            c = 1
            for i in dlist:
                try:
                    print(" %i of %i" % (c,len(dlist)))
                    os.system("scancel %s" % i)
                    c += 1
                except ValueError:
                    pass
            os.system("rm TMP_qdel")
        else:
            drange = drange.split("-")
            b = int(drange[0])
            e = int(drange[1])+1 # inclusive
            for i in range(b,e):
                print("scancel %s" % i)
                os.system("scancel %i" % i)
            
        print("Done!")

    #
    # Take a list of Job IDS in the following format:
    #   12345
    #   12345.cmgr01
    #
    def qdel2(self,klist):
        
        klist = open(klist).readlines()
        
        for i in klist:
            i = i.strip().split(".")[0]
            print(i)
            os.system("scancel %s" % i)
            
        print("Done!")

    def check_running(self,user):
        '''
        Note: Not tested with the new slurm scheduler, only change made was qstat to squeue -l
        '''
        
        os.system("squeue -l -u %s" % user)
        inp = open("TMP_check_running_log")
        inl = inp.readlines()
        countQ = 0
        countR = 0
        countE = 0
        countH = 0
        flagStart = 0
        flagUser = 0
        if user != "":
            flagUser = 1
            
        UD = {} # {username:{"R":x,"Q":x,"E":x,"H":x}}
        for i in inl:
            if i[:3] == "---":
                flagStart = 1
            elif flagStart:
                i = i.split(" ")
                tmp = []
                for j in i:
                    if j != "":
                        tmp.append(j)
                U = tmp[2]
                S = tmp[4]
                if U not in UD:
                    UD[U] = {S:1}
                elif S not in UD[U]:
                    UD[U][S] = 1
                else:
                    UD[U][S] +=1
                    
                if tmp[4] == "R":
                    countR += 1
                elif tmp[4] == "Q":
                    countQ += 1
                elif tmp[4] == "E":
                    countE += 1
                elif tmp[4] == "H":
                    countH += 1
                else:
                    print("UNK:",tmp[4])

        print("Running :",countR)
        print("Queueing:",countQ)
        print("Held    :",countH)
        print("Err     :",countE)
        flag = ['R','Q','H','E']
        if flagUser:
            if user in UD:
                counts = []
                for i in flag:
                    if i in UD[user]:
                        counts.append(UD[user][i])
                    else:
                        counts.append(0)
                print("User: %s,  %i running, %i in queue, %i held, %i with error" % \
                                (user,counts[0],counts[1],counts[2],counts[3]))
            else:
                print("USer: %s, no job..." % user)
        else:
            print("User: [R,Q,H,E,T]")
            users = UD.keys()
            users.sort()
            for i in users:
                print('%s: [' % i,)
                T = 0
                for j in flag:
                    if j in UD[i]:
                        print(UD[i][j],)
                        T += UD[i][j]
                    else:
                        print(0,)
                print('%i]' % T)

        os.system("rm TMP_check_running_log")
    
    def check_err(self,jobname):
        '''
        Note: Not updated for slurm
        '''
        files = os.listdir("./")
        
        # err log files
        efiles = []
        for i in files:
            if i.find(jobname) != -1 and i.find(".sh.e") != -1:
                efiles.append(i)
        efiles.sort()

        # go through err log files
        print("With error:")
        eidx = []
        for i in efiles:
            s = os.path.getsize("./%s" % i)
            try:
                idx = int(i[i.find(jobname)+len(jobname):i.find(".sh.e")])
                if s != 0:
                    inp = open(i)
                    inl = inp.readlines()
                    print(" %s:" % i,inl)
                    eidx.append(idx)
            except ValueError:
                print(" %s:" % i,"job index err")
                continue
        
        print("###############")
        print("Compile new command line. Note that this rely on the presence")
        print("of the shell script file in the same folder as the err log.")
        print("In addition, it just look for the last non-empty lines in .sh.")
        print("###############")
        print(" %i jobs failed" % len(eidx))
        oup = open("cmd_%s_witherr" % (jobname),"w")
        for i in eidx:
            # read the shell script
            inp = open("%s%i.sh" % (jobname,i))
            inl = inp.readlines()
            
            # Look for the last non-empty line
            inl.reverse()
            for j in inl:
                j = j.strip()
                if j != "":
                    oup.write("%s\n" % j)
                    break
            
        print("Command lines of jobs with err: cmd_%s_witherr" % jobname)

        print("Done!") 

    def rmlb(self,astr):
        astr = astr.strip()
        return astr

    def help(self):
        print("\nFunctions (-f):")
        print("    submit - create shell script and submit job bsaed on a file")
        print("       where each line is one job. NEED: c, OPT: w,m,J,p,ei,wd,mo,")
        print("       pre,post")
        print("    queue - submit jobs sequentially. NEED:c,u, OPT:s,n,r,w,m,")
        print("       J,p,e,k,wd,mo")
        print("    qdel  - delete jobs, NEED: r or u (if want to kill all)")
        print("    qdel2 - delete jobs based on a kill list (k)")
        print("    check_running - check how many job have R stat, OPT: u")
        print("    check_err - check any job error and compile a new cmd line")
        print("       file,  NEED: J")
        print("")
        print("Parameters:")
        print("    c - the command line file, one line per job, cmd start with #")
        print("        will not run.")
        print("    e - email, default ''. If given, FAIL, BEGIN, & END will be sent")
        print("    s - sec between qstat check, default 10")
        print("    nnode - number of nodes to request, default 1")
        print("    p - number of CPUs to use, default 1, same as -ncpu")
        print("    ngpu - number of GPUs to use, default 0")
        print("    n - number of jobs to submit at a time , default 10")
        print("    u - which user to monitor")
        print("    r - jobnum1-jobnum2, jobnum-, or all (default)")
        print("    w - wtime, in minutes. Defaulit 10 min.")
        print("    wd- working dir")
        print("    m - mem in GB, default 1")
        print("    mo- a file or a string with a list of modules to load, if a str,")
        print("        use ',' as delimiter")
        print("    J - job name")
        print("    pre - a file with additional commands common among all jobs to")
        print("        appended BEFORE the variable command line")
        print("    post- a file with additional commands common among all jobs to")
        print("        appended AFTER the variable command line")
        print("    o - keep all log files (1) or rename to tmp.o/tmp.e (0)")
        print("    k - a list of job ids to kill")
        print("    A - name of buy-in node")
        print("    array - Range if running an array job (i.e. 1-10)")
        print("    devnode - Specify which nodes to submit to (i.e.: intel16|intel18)")
        print("    i - Run in interactive mode, default = '' s")
        print("")
        sys.exit(0)

if __name__ == '__main__':

    qsub = qsub_hpc()
    f = c = u = e = wd = k = mo = pre = post = a = array = inta = devnode = ""
    s = n = 10
    w = 600
    m = o = p = nnode = 1
    ngpu = 0
    r = "*"
    J = "job"
        
    for i in range(1,len(sys.argv),2):
        if sys.argv[i] == "-f":
            f  = sys.argv[i+1]
        elif sys.argv[i] == "-c":
            c  = sys.argv[i+1]
        elif sys.argv[i] == "-e":
            e  = sys.argv[i+1]
        elif sys.argv[i] == "-s":
            s  = int(sys.argv[i+1])
        elif sys.argv[i] == "-n":
            n  = int(sys.argv[i+1])
        elif sys.argv[i] == "-nnode":
            nnode  = int(sys.argv[i+1])
        elif sys.argv[i] == "-p" or sys.argv[i] == "-ncpu":
            p  = int(sys.argv[i+1])
        elif sys.argv[i] == "-ngpu":
            ngpu  = int(sys.argv[i+1])
        elif sys.argv[i] == "-u":
            u  = sys.argv[i+1]
        elif sys.argv[i] == "-r":
            r  = sys.argv[i+1]
        elif sys.argv[i] == "-w":
            w  = int(sys.argv[i+1])*60
        elif sys.argv[i] == "-m":
            m  = int(sys.argv[i+1])
        elif sys.argv[i] == "-J":
            J  = sys.argv[i+1]
        elif sys.argv[i] == "-o":
            o  = int(sys.argv[i+1])
        elif sys.argv[i] == "-wd":
            wd = sys.argv[i+1]
        elif sys.argv[i] == "-k":
            k = sys.argv[i+1]
        elif sys.argv[i] == "-mo":
            mo = sys.argv[i+1]
        elif sys.argv[i] == "-pre":
            pre = sys.argv[i+1]
        elif sys.argv[i] == "-post":
            post = sys.argv[i+1]
        elif sys.argv[i] == "-A":
            a = sys.argv[i+1]
        elif sys.argv[i] == "-i":
            inta = sys.argv[i+1]
        elif sys.argv[i] == "-array":
            array = sys.argv[i+1]
        elif sys.argv[i] == "-devnode":
            devnode = sys.argv[i+1]
        else:
            print("UNKNOWN FLAG:",sys.argv[i])
            print("Do -h to get help.")
            sys.exit(0)

    if f == "submit":
        if "" in [c]:
            print("\nNeed cmd line file, user name\n")
            qsub.help()
        qsub.submit(c,0,w,m,J,p,e,o,a,inta,nnode,ngpu,array,devnode,wd,mo,pre,post)
    elif f == "queue":
        if "" in [c,u]:
            print("\nNeed cmd line file, user name\n")
            qsub.help()
        qsub.queue(c,s,n,u,r,w,m,J,p,e,o,a,inta,nnode,ngpu,array,devnode,wd,mo,pre,post)

    elif f == "qdel":
        if r == "" and u == "":
            print("\nNeed range or user name\n")
            qsub.help()
        qsub.qdel(u,r)
    elif f == "qdel2":
        if k == "":
            print("\nNeed kill list\n")
            qsub.help()
        qsub.qdel2(k)
    elif f == "check_running":
        qsub.check_running(u)
    elif f == "check_err":
        if J == "job":
            print("\nYou are using the default job name, make sure this is really what")
            print("you want to check... Go ahead anyway")
        qsub.check_err(J)
    else:
        print("\nUnknown function...\n")
        qsub.help()
