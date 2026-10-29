# Level ごとの domain decomposition と通信の隠蔽

作成: 2026-10-21 (Claude Code との作業メモ)

## 1. 現状の仕組み(変更前)

- `cpu_map(icell)` はセル単位の属性で、「そのセルの子 oct (grid) を所有する CPU」を表す。
  level `l` の grid の所有者は `cpu_map(father(igrid))`。
- `cpu_map` の値はセル中心の Hilbert key を `bound_key(0:ndomain)` の区間と比較して決める
  (`cmp_cpumap` / `cmp_ordering`)。`bound_key` は 1 組しかないので、
  **全 level で同じ区間分割**になる。
- `authorize_fine(ilevel)`: 他 CPU 所有の仮想セルは、その key 範囲
  `[order_min, order_max)` が自 CPU の区間と交われば「子を作ってよい」(flag2=1)。
  さらに `nexpand_bound` セル分だけ膨張させる。これが ghost (virtual) 領域の定義。
- 負荷分散 (`cmp_new_cpu_map`) は leaf cell のコスト (80 + 粒子数) × `niter_cost`
  (subcycle 回数) を Hilbert 順に並べて累積コストを `ndomain` 等分し、
  `bound_key2` を決める。つまり**全 level 合計の仕事量**だけを均す。

問題: level ごとの step は `MPI_ALLREDUCE` 等で同期するため、各 level の仕事量が
CPU 間で偏っていると、合計が均されていても待ち時間が出る。

## 2. 変更方針 (Phase 1: level ごとの分割)

### データ構造

- `bound_key(0:ndomain)` → `bound_key(0:ndomain, 1:nlevelmax)`。`bound_key2` も同様。
  第 2 添字は **grid の level** (= 親セルの level + 1)。
  - coarse cell (level 0, `1:ncoarse`) の `cpu_map` は `bound_key(:,1)` で決める。
  - level `l` のセルの `cpu_map` は `bound_key(:,l+1)` で決める。
- `cmp_cpumap(x,c,nn)` → `cmp_cpumap(x,c,nn,ilevel)` (`ilevel` = 作られる grid の level)。
  `cmp_dommap` も同様。
- 点 (位置) から「その点を含む leaf cell の所有者」を求める用途
  (`sink_particle.f90`, `flag_formation_sites.f90`) には新ルーチン
  `cmp_cpumap_leaf(x,c,nn)` を用意し、ローカル木を辿って leaf を見つけ、
  その grid の所有者を返す。自分の leaf なら `myid`、そうでなければ `myid` 以外。
- 粒子の初期配置 (`init_part`, `init_tracer`) は level 1 の grid に付けるので
  `cmp_cpumap(...,1)`。

### authorize_fine

level `l` の仮想セルは、`m = l+1 .. nlevelmax` のいずれかの `bound_key(:,m)` 区間
(自 CPU 分) と key 範囲が交われば authorize する。
理由: 自 CPU が level `m` の grid を所有するなら、その祖先 (level l のセル、
level l-1 のセル、…) はすべてローカルに存在しなければならないため
(`father(igrid)` はローカル index)。Hilbert key の範囲は階層的に入れ子なので、
この判定で祖先の存在が帰納的に保証される。`balance=.true.` 時は `bound_key2` でも同じ判定。

これが「通信が増える」主因: level `l` の ghost 領域が、自 CPU の
level `m>l` の担当領域の footprint をすべて含むようになる。
各 level の分割が似た区間になれば増加は小さく、大きくずれれば増える。
levelmin 未満の level は levelmin と同じ分割にして、無駄な ghost を避ける。

### cmp_new_cpu_map

level `m` ごとに独立に:
- 自 CPU 所有の level `m` の grid それぞれについて、親セル中心の key と
  コスト `twotondim*(80 + numbp(igrid)) * niter_cost(m)` を集める。
- Hilbert 順にソートして累積コストを `ndomain` 等分 → `bound_key2(:,m)`。
- 全体で grid がゼロの level は 1 つ粗い level の分割をコピー。
- `balance_per_level=.false.` (namelist, 既定値) なら従来通り 1 組だけ計算して
  全 level にコピー → 旧挙動と同じ。

### 入出力

- `amr_XXXXX.outYYYYY`: 従来の `bound_key(0:ndomain)` レコードには **全 level の
  hull** (区間の min/max) を書く (後処理ツールの互換性のため)。
  ファイル末尾に `bound_key(0:ndomain,1:nlevelmax)` レコードを追加する。
- restart 読み込み: 末尾レコードが読めればそれを使う。読めなければ (旧形式)
  1 組の値を全 level にコピー。
- `info_XXXXX.txt` の DOMAIN 表も hull。level ごとの表は `info_levelkeys_XXXXX.txt`
  に別途書く。
- `ordering='bisection'` との併用は不可 (エラー停止)。

## 3. 通信の隠蔽 (Phase 2)

RAMSES の通信は `make_virtual_fine_dp` (owner → ghost) と
`make_virtual_reverse_dp` (ghost → owner, 加算) の 2 種類で、どちらも
pack → ISEND/IRECV → WAITALL → unpack のブロッキング呼び出し。変数ごとに 1 メッセージ。

### 3a. メッセージの集約

`make_virtual_fine_dp_nv(xx, nvtot, iv1, iv2, ilevel)` /
`make_virtual_reverse_dp_nv(...)` を追加し、`uold(:,1:nvar)`、`rtuold(:,1:nrtvar)`
などを 1 メッセージ/相手で送る。NGROUPS=7 の RT では 28 変数 → 1 メッセージ。
レイテンシ律速の部分がそのまま減る。

### 3b. reverse 通信の非同期化 (計算との重ね合わせ)

`godfine1` は自 CPU の active セルと、粗い level (`l-1`) の隣接セルの `unew` にしか
書かない。level `l` の仮想セルへの寄与は level `l+1` の step (flux 補正) で
既に書き込まれており、`godunov_fine(l)` の呼び出し前に確定している。したがって

```
call make_virtual_reverse_start(unew, ilevel)   ! pack + ISEND/IRECV
call godunov_fine(ilevel)
call make_virtual_reverse_end  (unew, ilevel)   ! WAITALL + 加算
```

とすれば reverse 通信が hydro 計算に完全に隠れる。加算の順序は従来と同じ
(受信分を godunov の後に加算) なので結果は bit 単位で同じ。
RT (`rt_godunov_fine` + `add_rt_sources` は active セルのみ更新) も同様。

### 3c. forward 通信 (uold の ghost 更新)

step 末の `make_virtual_fine_dp(uold)` は直後の `flag_fine` (勾配判定) で必要なため、
重ねられる計算が少ない。候補は cooling_fine を「送信対象 grid → それ以外」の順に
分ける方法だが、star_formation の順序も動くので Phase 2 では見送り、
集約 (3a) のみ適用する。

## 4. 検証 (2026-10-21, Mac 上 4 ランク, gfortran 15)

| テスト | 結果 |
|---|---|
| Sedov AMR (levelmin 5, levelmax 7, nremap 2), `balance_per_level=.false.` | hydro 出力が変更前バイナリと **bit 一致** (Phase 1+2 込み) |
| 同上, `balance_per_level=.true.` | 完走。各 level の grid 数が min/max ±1 に均等化。保存量・min/max は単一分割と丸め誤差内で一致 |
| 同上 + `poisson=.true.` (multigrid) | 単一分割は bit 一致、level ごとも完走 |
| Stromgren RT AMR (rt=.true., 7 群, cooling) | 単一分割で hydro/rt 出力が bit 一致 (非同期 reverse 通信でも順序不変)、level ごとも完走 |
| restart: level ごと出力から (flag on / off), 旧形式出力から (flag on / off) | 全て完走。flag off で level ごとファイルを読むと自動で on に切替 |
| LIGHT_MPI_COMM=1 ビルド (Stromgren RT) | 単一分割で bit 一致、level ごとも完走 |
| DICE 銀河 + RT (`tests/dice_rt_amr.nml`: 本番 `xd_rt.nml` を levelmax 8, 6 step にしたもの。星の放射, cloudy 金属冷却, nremap=1) | 単一分割で hydro/rt/part 出力が HEAD と **bit 一致**。level ごとも完走し level 6-8 が ±0.1% に均等化 (単一分割は ±3-4%)。ghost: level 6 4794→9185, level 7 5488→12253, level 8 3955→4713 (1 CPU 最大は level 7 で 1380→5411) |
| DICE 銀河 (dice_dw1.g2, 粒子 194 万, levelmin 6, levelmax 8, 本番構成 RT+DICE, SF on, neq_chem) | 単一分割で hydro/part/grav 出力が HEAD と **bit 一致**。level ごとも完走し、level 6-8 の grid 数が ±0.2% に均等化 (単一分割では偏る)。ただし未追跡の `patch/takashi/flag_utils.f90` (Mizuki さんの若い星 refinement) を含めると HEAD でも変更後でも segfault するので、このテストはそのファイルを上流版に差し替えて実施 |

### ghost (virtual) oct の増加

load balance 後に `Level l has N ghost octs (min, max, ghost/active)` を表示する。
Sedov テストでの比較 (level 5 は一様な base level, 6-7 は爆発波付近のみ):

| level | 単一分割 ghost (ghost/active) | level ごと ghost (ghost/active) |
|---|---|---|
| 5 | 2331 (0.57) | 3142 (0.77) |
| 6 | 440 (0.48) | 745 (0.81) |
| 7 | 528 (0.61) | 579 (0.67) |

Stromgren (光源が箱の角にあり細かい level が偏在) では level 5 の ghost が 1 CPU で
3073 (単一分割なら ~300) と大きく増えた。細かい level の担当領域が粗い level の
担当領域から離れるほど、粗い level の ghost が増える (設計上の trade-off)。
本番 (levelmax 13-15, 粗い level の oct 数は全体の数 %) では粗い level の ghost
増加はメモリ・通信量ともに小さいはずだが、実測で確認すること。

## 4b. 本番 (896 ランク, levelmax 15) での segfault と修正 (2026-10-27)

`balance_per_level=.true.` で restart した本番ジョブ (flat-328709) が 4 回目の coarse step の
途中、level 15 の `multigrid_fine` → `build_parent_comms_mg` (MPI_Alltoall を含む) で
SIGSEGV。単一分割 (flat-328710) は同じバイナリで問題なし。

原因: **セルの `cpu_map` をどの level の `bound_key` で決めるかが 2 箇所で食い違っていた。**

- 設計 (`cmp_new_cpu_map`, `authorize_fine`, `cmp_new_bound_key_per_level`):
  level `l` のセルの `cpu_map` = そのセルの子 oct (level `l+1`) の所有者 = `bound_key(:,l+1)`。
- `refine_utils.f90` の `make_grid_fine(..., ilevel, ...)` (`ilevel` は新しく作る grid の level)
  では新 grid のセルに `cmp_cpumap(xx,cc,nn,ilevel)` を使っていた → `bound_key(:,ilevel)`
  で 1 level ずれる (`make_grid_coarse` も同様に 1 → 2)。

単一分割では全 level の `bound_key` が同じなので影響なし (bit 一致していた理由)。
level ごとの分割では、load balance 直後は `cmp_new_cpu_map` が全セルを正しく計算するが、
その後 refinement で作られた grid のセルだけ所有者が 1 level ずれる。すると
「`cpu_map` 上の所有者 p' が、その親セルを authorize しておらず grid を作っていない」
ケースが生じ、他 CPU はその grid を p' の仮想 grid として持つのに p' 側に実体がない。
`build_comm` で p' が受け取る grid index が 0 になり、

- ghost 交換 (`make_virtual_fine_dp*`) は `xx(iskip+0)` を黙って送る (物理量が壊れる)、
- multigrid は `lookup_mg=0` の grid を owner に要求 → owner 側で `father(0)` 等を参照して segfault。

4-8 ランクでは level 間の区間のずれが小さく、たまたま踏まなかった。

修正:
- `refine_utils.f90`: `make_grid_fine` で `cmp_cpumap(xx,cc,nn,ilevel+1)`、
  `make_grid_coarse` で `cmp_cpumap(xx,cc,1,2)`。
- `virtual_boundaries.f90` `build_comm`: 受信した emission の grid index に 0 があれば
  メッセージを出して `MPI_ABORT` する sanity check を追加 (所有者と authorize の不整合を
  黙って通さない)。コストは ghost 数の整数ループ 1 回。

### 4b-2. 修正後も本番で不整合が残る (2026-10-28, 未解決)

修正版 + `balance_level_min=13, balance_tolerance=0` の本番ジョブ (flat-329902) は 4 coarse step を
正常に進めた後、5 回目の途中 (level 13 の step 終了直後 = `amr_step(13,2)` 冒頭の refinement ループ、
`refine_fine(13)` → `refine_fine(14)` → `build_comm(15)`) で sanity check が発火した:

    build_comm: cpu 622 level 15: 27 virtual grids on other cpus have no owner copy here

つまり「他 CPU が cpu 622 の仮想 grid として持つ level 15 の grid 27 個が cpu 622 に無い」。
27 = 3×3×3 は `smooth_fine` の 1 セル buffer (孤立した flag の周り) と同じ数で、level 14 の
セル 1 個の周囲が他 CPU では refine され 622 では refine されなかった (または 622 がその grid
を持っていない) ことを示唆する。手元 (8 / 32 ランク, levelmax 8, 60 step, バグあり版でも) では再現しない。

コードを読み直した範囲では矛盾が見つからない (確認した点):
- `make_grid_fine` / `cmp_new_cpu_map` / `authorize_fine` の key 規約は一致 (セル中心の key に
  `bound_key(:,l+1)`、authorize はセルの key 範囲と全ての細かい level の区間の交差)。
  `cmp_ordering` は 2^(nlevelmax+1) 解像度なので丸めで隣セルに落ちることはない。
- 仮想セルの `flag1` は `authorize_fine` 末尾で所有者から復元、`cpu_map` は refinement ループで
  `build_comm(i)` 直後に所有者から同期されるので、同じ grid を持つ CPU 間で refine の判断は揃うはず。
- `cmp_new_bound_key_per_level` の境界は各 level で単調 (level 13-15 の grid 数 min/max が
  896 ± 数個なので map 自体は正常)。
- 整数通信ルーチン (`make_virtual_fine_int` / `make_virtual_reverse_int`) は上流と同一。

そこで `build_comm` の sanity check を**集団診断**に拡張した (`debug_missing_owner`):
`nbad` を `MPI_ALLREDUCE` し (1 int / 呼び出し)、1 つでもあれば
- 所有者は受け取った index 列を送り元に送り返す (echo)。
- 送り元 (holder) は index が 0 の仮想 grid について `MISSING holder ...` 行を出す:
  grid の `xg`、親セルの `cpu_map`/`flag1`/`son`、親 grid の所有者、親セル中心の
  `cmp_cpumap`、所有者がその親セルを authorize するか (level ごと)。位置を所有者に送る。
- 所有者は `OWNER ...` 行で、その位置の木を level 0 から辿り、各 level の grid の所有者・セルの
  `cpu_map`/`flag1`/`son`、自分がそのセルを authorize するか (level ごと) を出す。
  どの level から先が無いか、あるいは親セルはあるのに `son=0` なのかが分かる。
- その後 `MPI_ABORT`。通常時のコストは allreduce 1 回だけ。LIGHT_MPI_COMM ビルドでは診断なしで abort。
- `-DTEST_MISSING_OWNER` を付けてビルドすると cpu 3 が level levelmin+2 の coarse step 2 で
  わざと index を 0 にして診断経路を通す (動作確認済, 8 ランク)。

次の本番ジョブ (同じ設定) でこの出力を見て原因を特定する。`flag2(work)` は `build_comm` の
作業領域なので authorize の値ではない (authorize の判定は別途計算した T/F を見る)。

**診断の初出力 (flat-330537, 2026-10-28)**: 古い namelist で意図せず per-level (`balance_level_min`,
`balance_tolerance` は既定) で走ったジョブ。3 coarse step (569 / 424 / 493 s) 正常の後、4 step 目
(7579) に `build_comm: cpu 660 level 15: 27 virtual grids ... no owner copy` で発火。同時に
cpu 499 についても holder 506 が 8 grids を報告。診断行の要点:

- holder は cpu 664 (18 grids) と 665 (9 grids)。位置は (1.5137, 1.4676, 1.5005) 付近の 3×3×3 ブロック
  (level 15 の grid 間隔 2^-14)。全て `cmp_cpumap(centre,15)=660` で、660 はその親セルを level 15 で
  authorize する (T)。
- **holder 側の親セルが不正**: 親セル (level 14 grid 18886 / 17201 等のセル) の `cpu_map=0`,
  `flag1=0`, `son=0`、親 grid の所有者 (`cpu_map(father(fgrid))`) が 5242 (> ncpu=896)。
  仮想 grid の `father` が指すセルは生きた grid のセルに見えない (son であるはずの grid があるのに
  `son=0`)。holder の reception list にある level 15 grid が、親 grid が kill された後の stale な
  grid か、`father(igrid)` が壊れている可能性がある。
- **owner 660 側**: level 13 までは木があり (level 13 grid 37433 は owner 672 の仮想コピー)、
  そのセル ind 1/3/5/7 (`cpu_map=664` = level 14 oct の所有者) で `son=0`、level 14 の grid が無い。
  660 は level 14 を持たないので level 15 を所有できない。authorize 判定は `F T F` と出ているが、
  ログの診断行の形式が HEAD と違う (L 値が 1 つ多く `cmp_cpumap` が無い) ので、転送時点のコードが
  古いと思われ、順序の解釈は確定できない。

再開するなら次の一手: (1) holder 側で `reception(icpu,15)%igrid` の grid が本当に生きているか
(その親 grid が active/reception のどれかに含まれるか) を診断に足す。(2) owner 側の `flag1` は
`authorize_fine` の作業領域で信頼できないので、`refine_fine(13)` の時点で 660 がそのセルを
flag/authorize したかを別配列に記録する。per-level は 4c の結論で見送り中なので、ここまでで保留。

## 4c. 本番ログから見た性能 (2026-10-27) と ghost 抑制の設計

**結論 (2026-10-28): この問題・このランク数では per-level 分割は時間でもメモリでも得をしないので見送り**
(詳細は節末の表と説明)。


flat-328710 (単一分割) と flat-328709 (level ごと, バグあり) の比較。ghost は load balance 直後の値。

| | 単一分割 | level ごと |
|---|---|---|
| ghost oct 合計 (active 合計 521 万) | 298 万 (0.57) | 2076 万 (3.98) |
| level 15 / 14 の ghost/active | 0.32 / 0.53 | 1.51 / 3.19 |
| level 6-11 の ghost/active | 0.7-1.9 | 13.6-17.2 |
| grid メモリ | 27% | 39% |
| coarse step あたり | 307-432 s | 424-492 s |

**ghost が増える仕組み** (単一分割では起きない項):

1. 自分の level `l` grid のセルは全部 authorize される (upload に必要)。level `l` の担当領域 A_l と
   level `l+1` の担当領域 B_{l+1} が別の場所だと、A_l にある level `l+1` の oct を全部 ghost として
   持つ。この項だけで ghost/active ≈ 1 (細かい level では N_{l+1}/N_l × 自分の level-l oct 数)。
2. 粗い level では細かい level の担当領域の footprint (祖先) に 3×3×3 の膨張がかかる。
   896 CPU だと level 6 の自分の oct は 37 個しかないのに footprint の膨張で ~600 個の ghost。

**均等化の利益の見積り** (single ログの `Level l has N grids (min,max,avg)` から, 1 CPU の
最大 oct 数 × step 数の総和 = 計算のクリティカルパス):
単一分割 3.2M oct-step, 完全均等なら 1.8M (level 15 の max/avg 1.39, 14: 1.76, 13: 4.3)。
つまり計算は最大 ~45% 短縮できるが、ghost が 7 倍になると通信がそれを食い潰す
(観測は 1.2-1.4 倍遅い)。ghost を 2 倍以内に抑えつつ均等性の大部分を保つ必要がある。

**対策 (実装済)**: 2 つの namelist パラメータ (`&RUN_PARAMS`)。

- `balance_level_min` (既定 0 = levelmin): この level より粗い level は独立に均等化せず、
  この level の分割をコピーする。粗い level の footprint 膨張 (項 2) を消す。
- `balance_tolerance` τ (既定 0): 最細 level は厳密に均等化し、そこから粗い方へ順に、各境界を
  「その level の累積コストが目標 i·T/ndomain の ±(τ/2)·T/ndomain に収まる oct の範囲」の中で
  1 つ細かい level の境界に最も近い位置に置く (範囲内に細かい level の境界があればそれと一致
  させる)。境界が一致した所では項 1 が消える。両側の境界が逆向きにずれても domain のコストは
  平均の 1±τ に収まる。τ ≤ 0.9 にクランプ (窓が重ならない条件)。
  単一分割の level 13-14 の max/avg (4.3, 1.76) に比べ τ=0.5 でも均等。

実装は `load_balance.f90` の `cmp_new_bound_key_per_level` 後半。各 CPU がソート済みの
自分の oct を走査し、境界 j の窓に入る最初/最後の oct の key (klo/khi) と厳密境界を求め、
`MPI_ALLREDUCE` (MIN/MAX/SUM) で集約してから、細かい level の境界と比較して決める。
空の level は細かい level の分割をコピー、最細 level より上は最細 level をコピー。

**手元の検証 (2026-10-27, DICE RT テスト levelmin 6 / levelmax 8, 8 ランク, 12 step, bounds check 付き)**。
最後の load balance 後の ghost oct 数と、1 ランクあたりの grid 数 (min / max, 平均は level 6: 4096, 7: 4500, 8: 2854):

| 設定 | ghost L6 / L7 / L8 | grids L6 min/max | grids L7 min/max | grids L8 min/max |
|---|---|---|---|---|
| 単一分割 | 8112 / 8507 / 6214 | 3808 / 4294 | 4448 / 4539 | 2735 / 3051 |
| per-level τ=0 | 27284 / 25379 / 8160 | 4087 / 4102 | 4488 / 4510 | 2849 / 2864 |
| τ=0.3 | 19379 / 23297 / 8128 | 3328 / 4829 | 3704 / 5294 | 2849 / 2864 |
| τ=0.6 | 15488 / 20455 / 8125 | 2860 / 5578 | 2917 / 6089 | 2849 / 2864 |
| balance_level_min=8 | 10146 / 14266 / 7932 | 1468 / 7429 | 2533 / 6756 | 2849 / 2864 |
| balance_level_min=7, τ=0.3 | 11187 / 21048 / 8128 | 1573 / 7187 | 3704 / 5294 | 2849 / 2864 |

τ の不均衡の上限 (1±τ) は守られ、ghost は τ とともに単調に減る。この小さなテストでは単一分割が
もともとほぼ均等なので trade-off は見えにくい (本番では単一分割の level 13/14 が 4.3/1.76 倍)。
τ=0 は修正直後の per-level 出力と bit 一致、単一分割は変更前と bit 一致 (hydro/part/grav/rt/amr)。

**使い方の目安**: まず `balance_level_min=13, balance_tolerance=0.3` 程度から始め、load balance 後の
`Level l has N ghost octs` と `Level l has N grids (min,max,avg)` を単一分割と比べる。
ghost 合計が単一分割の 2 倍以内で level 13-15 の max/avg が 1.3 以下なら勝ち目がある。

**本番での最初の結果 (2026-10-28, 896 ランク, 同じ restart)**。4 本のジョブ:
329902 (`balance_level_min=13, τ=0`), 329906 (13, τ=0.3), 329907 (13, τ=0.6), 329908 (14, τ=0)。

| | 単一分割 | 13, τ=0 | 13, τ=0.3 | 13, τ=0.6 | 14, τ=0 |
|---|---|---|---|---|---|
| ghost oct 合計 (万) | 298 | 1490 | 1478 | 1442 | 1023 |
| 内 level 13/14/15 (M) | 0.84/0.90/0.66 | 3.2/5.5/3.1 | 3.2/5.4/3.1 | 3.1/5.3/3.0 | 1.9/3.3/3.1 |
| level 13/14/15 の max/avg | 4.25/1.81/1.38 | 1.00/1.00/1.00 | 1.30/1.30/1.00 | 1.49/1.60/1.00 | 3.74/1.00/1.00 |
| grid メモリ | 27% | 38% | 39% | 38% | 36% |
| 結末 | 正常 | 5 step 目で 4b-2 の abort | 4 step 目で同じ abort (cpu 622, 27 grids) | 3 step 正常, 4 step 目で時間切れ | (3 step 目途中までのログ) |

**coarse step の時間は比較できない**。step 時間は level 15 の RT subcycle 回数でほぼ決まり
(全ジョブで time ≈ 165-190 s + 0.39 s × [level 15 の RT subcycle 総数]; 1 coarse step = 256 fine step)、
その回数は dt に依存し、dt は星形成 (CPU ごとの乱数) で物理が分岐するためジョブごとに違う
(例: 3 step 目 単一分割 364 回 307 s, τ=0 256 回 271 s, τ=0.6 664 回 444 s)。
RT 以外の部分 (切片) は単一分割 165-177 s, per-level 170-188 s で差はなく、
**RT subcycle 1 回あたりの 0.39 s も分割によらない**: level 15 を完全に均等化しても RT の時間が
変わらないので、896 ランクでの RT step は CPU あたりの仕事量ではなく通信・同期 (allreduce,
ghost 交換の latency) で決まっている可能性が高い。均等化で稼げる余地は切片の一部だけ。

ghost については `balance_level_min` が効き (粗い level 890 万 → 310 万, 14 なら level 13 も消えて
1020 万)、τ は level 13-14 の ghost をほとんど減らさない (τ=0.6 でも 13-15 で 3%)。項 1 (自分の
grid の子が他人の grid) は境界が完全に一致しない限り残るので、τ の効果は小さい。

## 5. 使い方

```
&RUN_PARAMS
balance_per_level=.true.   ! 既定値 .false. (従来と完全に同じ挙動)
balance_level_min=13       ! 既定 0 (=levelmin): これより粗い level はこの level の分割をコピー
balance_tolerance=0.3      ! 既定 0: level ごとの許容不均衡。境界を細かい level に寄せて ghost を減らす (4c 節)
/
```

- `ordering='bisection'` とは併用不可。
- 出力: `info_levelkeys_XXXXX.txt` に level ごとの Hilbert key 区間。
  `info_XXXXX.txt` の DOMAIN 表は level ごと分割時には全域 [min,max] を書く
  (後処理ツールが Hilbert key でファイルを選ぶ場合に取りこぼさないため)。
- `amr_XXXXX.outYYYYY` の末尾に 2 レコード (`nlevelmax`, `bound_key(0:ndomain,1:nlevelmax)`)
  が追加される。旧コードでこのファイルを restart することはできない
  (level 1 の分割だけが読まれ、cpu_map と矛盾する)。
- Phase 2 (集約 + 非同期 reverse) は常に有効。`make_virtual_fine_dp_nv`,
  `make_virtual_reverse_dp_nv/_start/_end` (virtual_boundaries.f90)。

## 6. 変更ファイル (すべて patch/takashi 以下)

- 新規コピー+変更: load_balance.f90, refine_utils.f90, sink_particle.f90,
  flag_formation_sites.f90, init_tracer.f90
- 変更: amr_parameters.f90, amr_commons.f90, read_params.f90, to_dice/read_params.f90,
  init_amr.f90, output_amr.f90, virtual_boundaries.f90, amr_step.f90,
  init_part.f90, to_dice/init_part.f90
- ついでの修正: flag_utils.f90 (Mizuki さんの若い星まわり refinement ブロック)。
  `son(nbors_father_cells(1,ind))` で隣接 father cell が無い (index 0) ときに
  `son(0)` を参照する箇所をガードした。ただし同じブロックで `ind_grid(i)` が
  1..27 の範囲外になるケース (粒子が grid の外にある?) も bounds check で検出
  されており、こちらは未修正。`star=.true.` で 44 Myr より若い星があると踏む。
