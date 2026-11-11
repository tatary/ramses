#!/usr/bin/ruby

=begin
Hybrid MPI + OpenMP submit script (Slurm, Intel oneAPI).
Derived from run_partial_intel.rb: the node allocation is the same
(NumThreadsPerProc cores per MPI task, 112 cores per node), but the extra
cores are now used by OpenMP threads instead of being left idle.
The binary must be built with OPENMP=1 (bin/Makefile.xd_intel).

Case (1): Restart = "No"                       -> submit with xd_rt.nml
Case (2): Restart = "Yes", RestartFromQueue="No"  -> submit with xd_rt_rs.nml
Case (3): Restart = "Yes", RestartFromQueue="Yes" -> same, with -d afterany:RestartID
=end

#Select your system.
System = "xd"

#Restart = "Yes"
Restart = "No"
RestartFromQueue = "No"
#RestartFromQueue = "Yes"
#RestartID = "319856"

NumberofQueue=1
BinaryName1 = "./ramses/bin/ramses3d"
BinaryName2 = "./ramses/bin/ramses3d"
NumNodes=32
CoresPerNode=112
CoresPerSocket=56
ThreadsPerCore=2                                 # hardware threads per core: Slurm counts --cpus-per-task in these
                                                 # (339635: --cpus-per-task=4 gave 2 cores + their HT siblings per task)
CoresPerTask=4                                   # physical cores reserved per MPI task
NumThreadsPerProc=4                              # OpenMP threads per MPI task (<= CoresPerTask; 1 = flat run with the same placement)
UseNoMultithread = true                          # true: --hint=nomultithread, Slurm then counts cores and leaves
                                                 #       the HT siblings unused (cpus-per-task = CoresPerTask)
                                                 # false: cpus-per-task = CoresPerTask*ThreadsPerCore (whole cores incl. HT)
CpusPerTask = UseNoMultithread ? CoresPerTask : CoresPerTask * ThreadsPerCore
NumProc = NumNodes * CoresPerNode / CoresPerTask
NumTasksPerNode = CoresPerNode / CoresPerTask
NumTasksPerSocket = CoresPerSocket / CoresPerTask
OmpStackSize = "64m"                             # godfine1+unsplit need ~8 MB per thread
KmpBlockTime = "0"                               # ms a worker spins after a parallel region (Intel default 200):
                                                 # 0 lets idle workers sleep so the serial parts (poisson...) run at
                                                 # single-core turbo; "200" restores the default (TODO.md 2d-3)

RunName = "rt_hybrid"
StartUp = "srun"
#QName="M-large-a"
QName="M-test-a"
#QName="P-test-a"
#QName="P-large-a"
PinDebug = true                                  # I_MPI_DEBUG=4 (pin table from rank 0) and --cpu-bind=verbose (one line per task)

def write_script(str, snapnum, new)

    outputfile = open(str,"w")
    outputfile.print "#!/bin/bash\n"
    outputfile.print '#SBATCH --job-name=', RunName, "\n"
    outputfile.print '#SBATCH --partition=', QName, "\n"
    outputfile.print '#SBATCH --nodes=', NumNodes.to_s, "\n"
    outputfile.print '#SBATCH --ntasks=', NumProc.to_s, "\n"
    outputfile.print '#SBATCH --ntasks-per-node=', NumTasksPerNode.to_s, "\n"
    outputfile.print '#SBATCH --ntasks-per-socket=', NumTasksPerSocket.to_s, "\n"
    outputfile.print '#SBATCH --cpus-per-task=', CpusPerTask.to_s, "\n"
    if UseNoMultithread
        outputfile.print "#SBATCH --hint=nomultithread\n"
    end
    outputfile.print "#SBATCH -o flat-%J.out\n"
    #outputfile.print "#SBATCH --mem=461G\n"    # memory per node
    outputfile.print "#SBATCH --mem=115G\n"    # memory per node
    #outputfile.print "#SBATCH --time=24:00:00\n"
    outputfile.print "#SBATCH --time=00:30:00\n"
    outputfile.print "#SBATCH --mail-type=BEGIN,END\n"
    outputfile.print "#SBATCH --mail-user=t.t.okamoto@gmail.com\n"
    outputfile.print ". /work/opt/local/bin/enable-oneapi.sh \n"
    outputfile.print "module swap PrgEnv-oneapi/2024.2 PrgEnv-oneapi/2025.3 \n"
    outputfile.print "export SLURM_MPI_TYPE=pmi2 \n"
    outputfile.print "export LD_LIBRARY_PATH=/lib64:${LD_LIBRARY_PATH} \n"
    #outputfile.print "export FI_OFI_RXM_USE_SRX=0 \n"
    #outputfile.print "export FI_VERBS_PREFER_XRC=0 \n"
    outputfile.print "\n"
    outputfile.print "# --- hybrid MPI + OpenMP ---\n"
    # --cpus-per-task in the header is not inherited by srun on Slurm >= 22.05
    outputfile.print "export SRUN_CPUS_PER_TASK=", CpusPerTask.to_s, "\n"
    outputfile.print "export OMP_NUM_THREADS=", NumThreadsPerProc.to_s, "\n"
    outputfile.print "export MKL_NUM_THREADS=1\n"
    # worker-thread stacks: the OpenMP work arrays live on the stack (ulimit -s only
    # covers the main thread; with 'unlimited' glibc even falls back to 2 MB threads)
    outputfile.print "export OMP_STACKSIZE=", OmpStackSize, "\n"
    # keep each thread on its own core inside the task's cpu set (a place is a core
    # with both hardware threads, so NumThreadsPerProc <= CoresPerTask means no HT sharing)
    outputfile.print "export OMP_PLACES=cores\n"
    outputfile.print "export OMP_PROC_BIND=close\n"
    outputfile.print "export KMP_BLOCKTIME=", KmpBlockTime, "\n"
    # let Slurm pin each task to its CoresPerTask cores (--cpu-bind=cores below);
    # Intel MPI then stays inside that cpu set, so I_MPI_PIN_DOMAIN is not needed
    if PinDebug
        outputfile.print "export I_MPI_DEBUG=4\n"
    end
    outputfile.print "\n"
    outputfile.print "ulimit -s unlimited\n"
    outputfile.print "ulimit -c 0\n"
    outputfile.print "\n"
    outputfile.print "cd ${SLURM_SUBMIT_DIR}\n"
    cpubind = PinDebug ? "verbose,cores" : "cores"
    # --hint is repeated on srun: like --cpus-per-task it is not reliably inherited from the header
    hint = UseNoMultithread ? " --hint=nomultithread" : ""
    srun = StartUp + " --cpus-per-task=" + CpusPerTask.to_s + hint + " --cpu-bind=" + cpubind
    if new == 1
        outputfile.print srun, " ", BinaryName1, " xd_rt.nml \n"
    else
        outputfile.print srun, " ", BinaryName2, " xd_rt_rs.nml \n"
    end

    outputfile.flush
    outputfile.close
end

def Attention
  if Restart == "Yes"
    if RestartFromQueue == "Yes"
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
    str = result.chomp!
    qid = str[/\d+$/]
    write_script("tmp_", i,2)
    `chmod 744 ./tmp_`
    command = "sbatch -d afterany:" + qid.to_s + " ./tmp_"
    result = `#{command}`
    print "Submit ",i+1,"-th script:",result.chomp,"\n"
    `rm -rf ./tmp_`
  end
}
