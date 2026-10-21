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

## 5. 使い方

```
&RUN_PARAMS
balance_per_level=.true.   ! 既定値 .false. (従来と完全に同じ挙動)
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
