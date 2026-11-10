#!/bin/bash
# OpenMP ビルドの実行と参照ランとの bit 比較 (引き継ぎ用の汎用版)。
# usage: run_omp_cmp.sh <binary> <ref_dir> <tag> <threads> <nranks> [namelist]
#   <ref_dir>  : 参照ラン (run.nml と output_0000?/ がある) のディレクトリ。
#   <tag>      : ref_dir の隣に作る実行ディレクトリ名。
#   namelist 省略時は ref_dir/run.nml をコピーする。
# 実行ディレクトリの親に dice/ UVB/ cloudy/ SED/ のシンボリックリンクが必要 (CLAUDE.md テスト節)。
# 比較は BITREPRO=1 でビルドしたバイナリ同士で行うこと (FMA 契約の差を避ける; TODO.md 2d-2)。
set -u
BIN=$1; REF=$(cd "$2" && pwd); TAG=$3; NT=$4; NP=$5; NML=${6:-$REF/run.nml}
cd "$(dirname "$REF")" && rm -rf "$TAG" && mkdir "$TAG" && cp "$NML" "$TAG/run.nml" && cd "$TAG"
ulimit -s unlimited 2>/dev/null || ulimit -s "$(ulimit -Hs)"
export OMP_NUM_THREADS=$NT OMP_STACKSIZE=128m
/usr/bin/time -p mpirun -np "$NP" -x OMP_NUM_THREADS -x OMP_STACKSIZE "$BIN" run.nml > log.txt 2> time.txt
echo "$TAG exit $? ranks $NP threads $NT $(grep real time.txt)"
n=0; d=0
for f in "$REF"/output_0000?/{amr,hydro,part,grav,rt}_*.out0000?; do
  [ -e "$f" ] || continue
  g=${f#$REF/}; n=$((n+1)); cmp -s "$f" "$g" || { d=$((d+1)); echo "DIFF $g"; }
done
echo "$TAG vs $(basename "$REF"): compared $n files, $d differ"
grep "hydro - godunov\|rt - godunov\|cooling  \|TOTAL" log.txt | tail -4 | cut -c1-100
