#!/usr/bin/ruby

=begin
This is a sequential queue submit script.
=end

#Select your system.
System = "xd"
#System = "opt"

#If Restart is "No", this script submit a job with "subfind.param" 
#If Restart is "Yes", this script submit a job with "subfind2.param"
#Case (1): RestartFromQueue = "No", simply submit a job
#Case (2): RestartFromQueue = "Yes", submit a job with -W depend=afterany:RestartID.  
Restart = "Yes"
#Restart = "No"
RestartFromQueue = "No"
#RestartFromQueue = "Yes"
#RestartID = "5943"

NumberofQueue=1
BinaryName1 = "./ramses_dice/bin/ramses3d"
BinaryName2 = "./ramses/bin/ramses3d"
NumNodes=8
NumThredsPerProc=1
NumProc = NumNodes * 112 / NumThredsPerProc
NumTasksPerNode=112/NumThredsPerProc
NumTasksPerSocket = 56 / NumThredsPerProc

#RunName = "rt_sn"
RunName = "test"
StartUp = "srun --hint=nomultithread"     # repeated on srun: not reliably inherited from the header
#QName="M-large-b"
#QName="M-test-b"
#QName="P-test-b"
QName="P-large-bp"
Lcommand='nodes=' + NumNodes.to_s

def write_script(str, snapnum, new)

    outputfile = open(str,"w")
    outputfile.print "#!/bin/bash\n"     
    outputfile.print '#SBATCH --job-name=', RunName, "\n"
    outputfile.print '#SBATCH --partition=', QName, "\n"
    outputfile.print '#SBATCH --nodes=', NumNodes.to_s, "\n"
    outputfile.print '#SBATCH --ntasks=', NumProc.to_s, "\n"
    outputfile.print '#SBATCH --ntasks-per-node=', NumTasksPerNode.to_s, "\n"
    outputfile.print '#SBATCH --ntasks-per-socket=', NumTasksPerSocket.to_s, "\n"
    outputfile.print "#SBATCH --hint=nomultithread\n"     # one task per core, HT siblings left unused
    outputfile.print "#SBATCH -o flat-%J.out\n" 
    outputfile.print "#SBATCH --mem=461G\n"    # ノード当たりのメモリ容量
    #outputfile.print "#SBATCH --mem=115G\n"    # ノード当たりのメモリ容量
    outputfile.print "#SBATCH --time=24:00:00\n"   
    outputfile.print "#SBATCH --mail-user=t.t.okamoto@gmail.com\n"
    #outputfile.print ". /work/opt/local/bin/enable-oneapi.sh \n"
    outputfile.print "source /work/opt/local/bin/enable-cpe.sh \n"
    outputfile.print "module switch PrgEnv-cray PrgEnv-gnu \n"
    outputfile.print "export SLURM_MPI_TYPE=pmi2 \n"
    outputfile.print "export LD_LIBRARY_PATH=/lib64:${LD_LIBRARY_PATH} \n"
    outputfile.print "export FI_OFI_RXM_USE_SRX=0 \n" 
    outputfile.print "export FI_VERBS_PREFER_XRC=0 \n"
    # a binary built with OPENMP=1 would otherwise start one thread per logical CPU in every rank
    outputfile.print "export OMP_NUM_THREADS=", NumThredsPerProc.to_s, "\n"
    outputfile.print "export OMP_STACKSIZE=64m\n"
    outputfile.print "ulimit -s unlimited\n"
    outputfile.print "\n"
    outputfile.print "cd ${SLURM_SUBMIT_DIR}\n"
    if new == 1
        outputfile.print StartUp, " ", BinaryName2, " " 
        outputfile.print "xd_rt.nml > log.", snapnum.to_s, " 2>&1" "\n"
    else
        outputfile.print StartUp, " ", BinaryName2, " " 
        outputfile.print "xd_rt_rs.nml > log.", snapnum.to_s, " 2>&1" "\n"
    end

    outputfile.flush
    outputfile.close
end

def Attention 
  if Restart == "Yes" 
    if RestartFromQueue == "YES" 
      if RestartID == "" 
        print "Incorrect RestartID\n" 
        exit(1) 
      end 
    end
  end
end

Attention()

result = 0
(NumberofQueue).times{ |i| 
  if i == 0 
    if Restart == "No" 
      write_script("tmp", i, 1) 
      `chmod 744 ./tmp` 
      result = `sbatch ./tmp` 
      print result
      print "Submit first script:",result.chomp,"\n" 
      `rm -rf ./tmp` 
    else 
      if RestartFromQueue == "No" 
        write_script("tmp_", i, 2) 
        `chmod 744 ./tmp_` 
        result = `sbatch ./tmp_` 
        print "Submit first script:",result.chomp,"\n" 
        `rm -rf ./tmp_` 
      else 
        write_script("tmp_", i, 2) 
        `chmod 744 ./tmp_` 
        command = "sbatch -d afterany:" + RestartID + " ./tmp_" 
        result = `#{command}` 
        print "Submit first script:",result.chomp,"\n" 
        `rm -rf ./tmp_` 
      end 
    end 
  else 
    print "hoge\n"
    print result, "\n" 
    print "moge\n"
    str = result.chomp!
    qid = str[/\d+$/]
    print qid, "\n"
    write_script("tmp_", i,2) 
    `chmod 744 ./tmp_` 
    command = "sbatch -d afterany:" + qid.to_s + " ./tmp_" 
    print "hoge\n"
    print command, "\n"
    print "moge\n"
    result = `#{command}` 
    print "Submit ",i+1,"-th script:",result.chomp,"\n" 
    `rm -rf ./tmp_` 
  end
}

