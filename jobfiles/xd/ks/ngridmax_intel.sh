#!/bin/bash
#SBATCH --job-name=ks_M10
#SBATCH --partition=P-large-a
#SBATCH --nodes=4
#SBATCH --hint=nomultithread
#SBATCH --ntasks=448
#SBATCH --ntasks-per-node=112
#SBATCH --ntasks-per-socket=56
#SBATCH --mem=461G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=t.t.okamoto@gmail.com
#
# patch/takashi 版の、Ramses-FG `jobfiles/xd/smbh/ngridmax_intel.sh` に対応する job script (2026-11-12)。
# 物理パラメータは ks の本番 (ngridmax_fire3.nml + sed 上書き) に揃えてある: namelist/ks/ngridmax_fire3.nml
# (takashi にまだ無い項目は !KS でコメント)。ビルドは bin/Makefile.xd_ks_intel (bin_ks/ で組む)。
# ks 側との違い: binary は ramses_ks3d、namelist は sed で上書きせずファイルをそのまま使い、
# nrestart と初期条件・テーブルのパスだけここで埋める。flat ラン (OpenMP なし; ks と同じ 448 ランク)。
# ハイブリッドで回すときは jobfiles/run_hybrid_intel.rb の設定 (--cpus-per-task, OMP_*) を足す。
. /work/opt/local/bin/enable-oneapi.sh
module swap PrgEnv-oneapi/2024.2 PrgEnv-oneapi/2025.3
module list
export SLURM_MPI_TYPE=pmi2
export LD_LIBRARY_PATH=/lib64:${LD_LIBRARY_PATH}
ulimit -a
ulimit -s unlimited
ulimit -c 0
set -x

NRESTART=-1                                      # restart data file (use latest data if NRESTART < 0)
RAMSESDIR=/home/$USER/ramses                     # patch/takashi のリポジトリ (ks 側は /home/$USER/Ramses-FG)
BINDIR=$RAMSESDIR/bin_ks                         # bin/Makefile.xd_ks_intel で組んだディレクトリ
RAMSES=ramses_ks3d
ICS=/work/$USER/smbh/ICs/4sigmaM10               # initial data directory (ks と同じ)
NMLBASE=namelist/ks/ngridmax_fire3.nml           # base namelist (ks の sed 上書き適用済み)
TABLES=/work/$USER/dwarf/ramses_tables           # SED / cloudy / UVB (xd_rt.nml と同じ場所)
SED=$TABLES/SED/bc03_Chabrier                    # stellar SED (ks と同じ bc03_Chabrier)
CHTAB=$TABLES/cloudy_cooling/cloudy_excess_metal_cooling.bin   # 方針 1: takashi の cloudy テーブル

rm -f ./$RAMSES && cp $BINDIR/$RAMSES ./ || { echo "ERROR: failed to install $RAMSES" >&2; exit 1; }
# namelist が相対パスで参照するもの (CHARACTER(80) に収めるため symlink)
ln -sfn $ICS ICS
ln -sfn $SED SED
mkdir -p cloudy && ln -sfn $CHTAB cloudy/cloudy_excess_metal_cooling.bin

if [ $NRESTART -lt 0 ]; then
    NRESTART=`ls -1d output* 2>/dev/null | wc -l`   # use latest data
    if [ $NRESTART -ge 1 ]; then
        NRESTART=-1                                 # use restart files
    else
        NRESTART=0
    fi
fi

NMLFILE=namelist.nml
cp $RAMSESDIR/$NMLBASE $NMLFILE
sed -i -e "s@^nrestart.*@nrestart=$NRESTART@g" $NMLFILE
sed -i -e "s@^walltime_hrs.*@walltime_hrs=23.5@g" $NMLFILE
sed -i -e "s@^restart_hrs.*@restart_hrs=1.0@g" $NMLFILE
sed -i -e "s@^initfile(1).*@initfile(1)='ICS/r4sigma.grafic_128'@g" $NMLFILE
sed -i -e "s@^initfile(2).*@initfile(2)='ICS/r4sigma.grafic_256'@g" $NMLFILE
sed -i -e "s@^initfile(3).*@initfile(3)='ICS/r4sigma.grafic_512'@g" $NMLFILE
sed -i -e "s@^initfile(4).*@initfile(4)='ICS/r4sigma.grafic_1024'@g" $NMLFILE
sed -i -e "s@^sed_dir.*@sed_dir='SED'@g" $NMLFILE
sed -i -e "s@^cloudy_metal_file.*@cloudy_metal_file='cloudy/cloudy_excess_metal_cooling.bin'@g" $NMLFILE
grep -q "^nrestart=$NRESTART" $NMLFILE || { echo "ERROR: nrestart not set in $NMLFILE" >&2; exit 1; }

NUM=`ls -1 stdout* 2>/dev/null | wc -l`
STDOUT="stdout"$((NUM))
date 2>&1 | tee -a $STDOUT
sync

# restart_files の退避 (ks と同じ: 前回の job が dump を書かなかったときは回さない)
newest_mtime() { find "$1" -type f -printf '%T@\n' 2>/dev/null | sort -n | tail -1 | cut -d. -f1; }
if [ -d restart_files ]; then
    live_t=$(newest_mtime restart_files)
    bak_t=$(newest_mtime restart_files.bak)
    if [ -d restart_files.bak ] && [ -n "$live_t" ] && [ -n "$bak_t" ] && [ "$live_t" -le "$bak_t" ]; then
        echo "restart_files has no files newer than restart_files.bak; skipping rotation (previous job probably wrote no dump)"
    else
        rm -rf restart_files.bak.tmp
        if [ -d restart_files.bak ]; then
            rm -rf restart_files.bak2
            mv restart_files.bak restart_files.bak2 || { echo "ERROR: failed to rotate restart_files.bak" >&2; exit 1; }
        fi
        cp -a restart_files restart_files.bak.tmp || { echo "ERROR: failed to copy restart_files (disk full?)" >&2; rm -rf restart_files.bak.tmp; exit 1; }
        mv restart_files.bak.tmp restart_files.bak || { echo "ERROR: failed to finalize restart_files.bak" >&2; exit 1; }
        echo "restart_files backed up to restart_files.bak (previous -> restart_files.bak2)"
    fi
fi

cd ${SLURM_SUBMIT_DIR}
srun --hint=nomultithread --cpu-bind=cores ./$RAMSES $NMLFILE 2>&1 | tee -a $STDOUT
date 2>&1 | tee -a $STDOUT
