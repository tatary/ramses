#!/bin/bash
#SBATCH --job-name=M10_music
#SBATCH --partition=P-large-a
#SBATCH --nodes=4
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-socket=8
#SBATCH --cpus-per-task=7
#SBATCH --hint=nomultithread
#SBATCH --mem=461G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=t.t.okamoto@gmail.com
#
# MUSIC 初期条件 music_M1e10_in_M5e11 (4sigmaM10 の後継) を patch/takashi の ks 設定で走らせる job script (2026-11-27)。
# hybrid: 4 ノード (元の 4sigmaM10 と同じ) x 16 ランク x 7 スレッド = 64 ランク (14 コアの NUMA ノードに 2 ランク)。
# 物理は ks の Pop II 比較ラン C_eps1 にそろえてある: namelist/ks/music_M1e10_in_M5e11.nml。
# ビルド (oneAPI, OpenMP):  mkdir -p bin_ks && cp bin/Makefile.xd_ks_intel bin_ks/Makefile && cd bin_ks && make -j1
# 実行ディレクトリで sbatch する。何本かまとめて投入してよい (run.lock/ で同時実行を防ぎ、2 本目以降は
# output_restart/ から再開する)。namelist は毎回リポジトリからコピーし、nrestart だけここで埋める。
. /work/opt/local/bin/enable-oneapi.sh
module swap PrgEnv-oneapi/2024.2 PrgEnv-oneapi/2025.3
module list
export SLURM_MPI_TYPE=pmi2
export LD_LIBRARY_PATH=/lib64:${LD_LIBRARY_PATH}

# --- hybrid MPI + OpenMP (jobfiles/xd/hybrid/run_hybrid_intel_nodep.rb と同じ設定) ---
NTHREADS=7
export SRUN_CPUS_PER_TASK=$NTHREADS     # --cpus-per-task は Slurm >= 22.05 では srun に引き継がれない
export OMP_NUM_THREADS=$NTHREADS
export MKL_NUM_THREADS=1
export OMP_STACKSIZE=128m               # OpenMP の作業配列は各スレッドの stack に乗る
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export KMP_BLOCKTIME=1

ulimit -a
ulimit -s unlimited
ulimit -c 0
set -x

RUN=music_M1e10_in_M5e11
RAMSESDIR=/home/$USER/ramses                     # patch/takashi のリポジトリ
BINDIR=$RAMSESDIR/bin_ks                         # bin/Makefile.xd_ks_intel で組んだディレクトリ
RAMSES=ramses_ks3d
NMLBASE=$RAMSESDIR/namelist/ks/$RUN.nml
ICS=/work/$USER/smbh/ICs/$RUN                    # MUSIC の出力 (level_008 ...; ics_ramses/ の下にあればそちらを使う)
TABLES=/work/$USER/dwarf/ramses_tables           # SED / cloudy (xd_rt.nml と同じ場所)
SED=$TABLES/SED/bc03_Chabrier                    # stellar SED (ks と同じ bc03_Chabrier)
CHTAB=$TABLES/cloudy_cooling/cloudy_excess_metal_cooling.bin   # takashi の cloudy テーブル (docs/merge_ks.md 方針 1)

cd ${SLURM_SUBMIT_DIR}

# 1 つのディレクトリで走る job は 1 本だけ: run.lock/jobid が持ち主。queue に居ない job の lock は引き継ぐ
if ! mkdir run.lock 2>/dev/null; then
    owner=$(cat run.lock/jobid 2>/dev/null)
    if [ -n "$owner" ] && squeue -h -j "$owner" -t RUNNING 2>/dev/null | grep -q .; then
        echo "job $owner is still running in this directory: exit"; exit 0
    fi
    echo "stale lock of job $owner: taking over"
fi
echo $SLURM_JOB_ID > run.lock/jobid
unlock() { if [ "$(cat run.lock/jobid 2>/dev/null)" = "$SLURM_JOB_ID" ]; then rm -rf run.lock; fi; }
die() { echo "ERROR: $*" >&2; unlock; exit 1; }

[ -d $ICS/ics_ramses/level_008 ] && ICS=$ICS/ics_ramses
[ -f $ICS/level_008/ic_deltab ] || die "no initial conditions in $ICS (level_008/ic_deltab)"
[ -d $SED ] || die "no SED directory $SED"
[ -f $CHTAB ] || die "no cloudy table $CHTAB"
[ -f $NMLBASE ] || die "no namelist $NMLBASE"

rm -f ./$RAMSES && cp $BINDIR/$RAMSES ./ || die "failed to install $RAMSES"
# namelist が相対パスで読むもの (CHARACTER(80) に収めるため symlink)
ln -sfn $ICS ICS
ln -sfn $SED SED
mkdir -p cloudy && ln -sfn $CHTAB cloudy/cloudy_excess_metal_cooling.bin

# restart: patch/takashi は restart_hrs ごとに output_restart/ を書き (前の組は output_restart.prev/ に残る)、
# nrestart=-1 でそこから再開する。
# 番号付きの output から再開するときは、投入前に実行ディレクトリで `echo 12 > restart_from` とする (output_00012 から)。
# restart_from はこの job が新しい output_restart/ を書くまで残すので、最初の dump の前に落ちても
# 次の job がまた output_00012 から始める (古い output_restart/ は使わない)。書けたら restart_from.done に改名する。
if [ -f restart_from ]; then
    NRESTART=$(tr -d ' \n' < restart_from)
    [[ "$NRESTART" =~ ^[0-9]+$ ]] && [ $((10#$NRESTART)) -ge 1 ] || die "restart_from must hold an output number (got '$NRESTART')"
    NRESTART=$((10#$NRESTART))
    RSDIR=$(printf "output_%05d" $NRESTART)
    [ -d $RSDIR ] || die "restart_from: no $RSDIR"
    echo "restart_from: restarting from $RSDIR"
elif [ -d output_restart ]; then NRESTART=-1; else NRESTART=0; fi

NMLFILE=namelist.nml
cp $NMLBASE $NMLFILE
sed -i -e "s@^nrestart.*@nrestart=$NRESTART@g" $NMLFILE
grep -q "^nrestart=$NRESTART\$" $NMLFILE || die "nrestart not set in $NMLFILE"

NUM=`ls -1 stdout* 2>/dev/null | wc -l`
STDOUT="stdout"$((NUM))
date 2>&1 | tee -a $STDOUT
sync

# --cpus-per-task と --hint は header から srun に確実には引き継がれないので繰り返す
srun --cpus-per-task=$NTHREADS --hint=nomultithread --cpu-bind=cores ./$RAMSES $NMLFILE 2>&1 | tee -a $STDOUT
date 2>&1 | tee -a $STDOUT
# restart_from は、この job が output_restart/ を書き終えたら役目を終える (次からは nrestart=-1)
if [ -f restart_from ] && [ -d output_restart ] && [ output_restart -nt restart_from ]; then
    mv -f restart_from restart_from.done
fi
unlock
