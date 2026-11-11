# やり残し・既知の問題 (引き継ぎメモ)

最終更新: 2026-10-29 (小野さんの flag_utils 修正の取り込み、RT ソルバの続き)

## 1. 既知のバグ

### 1a. (修正済 2026-10-29) `patch/takashi/flag_utils.f90` の若い星まわり refinement
- 2026-10-21 に `flag_utils.f90` → `flag_utils.f90.mizuki` に改名した。現在は上流の
  `amr/flag_utils.f90` がビルドされる。修正後に `.f90` に戻して `git add` する。
- 小野さんが `userflag_fine` に追加したブロック (`if(star)` … `agemax=0.044` Gyr)。
  44 Myr より若い星粒子を含む grid の、3×3×3 親セルの子 grid を強制的に細分化する。
- 問題 1: 隣接する親セルが存在しないとき `son(nbors_father_cells(1,ind)=0)` を参照
  → `son(0)` で segfault。**ガード済み** (`if(nbors_father_cells(1,ind_grid(i))>0)`)。
- 問題 2: `ind_grid(i)` (1..27 のはず) が範囲外になる。粒子位置から親セルの添字を
  求める計算 (`x0_part=xg-3*dx`, `dd_part=x_part+0.5`, `igg/igd`) が、粒子が grid の
  外に半セル以上出ている場合に 1..27 を超える。`move_fine` 直後 (粒子が grid を
  出た後、`make_tree_fine` 前) に `flag_fine` が呼ばれるのが原因と推測。**未修正**。
  bounds check 付きビルド (`-fbounds-check`) で DICE テスト (tests/dice_amr.nml) を
  流すと `flag_utils.f90` の `son(...)` 行で再現する。
- このファイルがあると HEAD でも変更後でも DICE 銀河が初回 load balance 直後に落ちる。
  develop (このファイルなし) は落ちない。修正するまでは外してビルドするのが安全。
- 修正したら `git add` して追跡下に置くこと。
- **2026-10-28 追記**: origin/mizuki (4112010f) の新版は問題 1 (`son(0)`) にガードが入り ghost の flag を
  `make_virtual_reverse_int` で戻す同期も加わったが、**問題 2 は未修正**。`tests/dice_amr.nml` を
  bounds check 付き 8 ランクで流すと最初の coarse step で `userflag_fine` の
  `nbors_father_cells(i,ind_grid_nbor(ix))` (533 行) が添字 28-33 (上限 27) になり、続く
  `flag1(ind_cell_nbor(...))` (540 行) が負の添字になる。これは takashi に移植した場合だけでなく、
  mizuki ブランチそのもの (`patch/mizuki`, SOLVER=mhd, `riemann='hlld'` に変更した同じ namelist) でも
  同じ行・同じ添字で再現する。原因は粒子が親セルの 3×3×3 の外にある場合 (`igg_loc/igd_loc`=3)。
  新版は `patch/takashi/flag_utils.f90.mizuki` として置いてある (ビルド対象外)。
- **2026-10-29 解決**: origin/mizuki acb88172 (小野さん) が `userflag_fine` を書き直した。ilevel の星と
  levelmax (ilevel+1) の星の走査を分け、後者を vector ループの外で nvector ずつ処理するようにした
  (粒子数は `numbp` で回す)。問題 2 の再現手順 (tests/dice_amr.nml, 8 ランク, `-fbounds-check`) で
  16 fine step 完走し、添字範囲外は出なくなった。そのまま `patch/takashi/flag_utils.f90` にコピーして
  ビルドに戻した (`.f90.mizuki` は削除)。上流 `amr/flag_utils.f90` との差は小野さんの追加のみ。
  これにより 44 Myr より若い星粒子のまわりの細分化が本番にも入る (物理結果が変わる)。

### 1b. `to_dice/init_part.f90`: `zp` の確保条件
- `zp` (粒子金属量) は `if(star.or.sink)` の中で `if(metal)` のときだけ確保されるが、
  `to_dice/init_flow_fine.f90` は `metal` だけ見て `zp` を読む。`eps_star=0` で
  `metal=.true.` だと未確保配列を読む (本番 namelist では起きない)。
  `if(metal)` の確保を `star.or.sink` の外に出せばよい。

### 1c. `rt_init.f90` ~316 行: `rt=.false.` かつ `neq_chem=.false.` で落ちる
- NIONS=4 ビルドで ion の添字 (`ixHII` 等) が設定されず `group_csn(…,0)` を参照。
  本番では常に rt=.true. なので実害なし。テスト時の罠。

### 1d. (コード削除 2026-10-29) balance_per_level の本番での所有者不整合
- per-level 分割のコードごと削除したので消滅。経緯は `docs/level_domain_decomposition.md` 4b/4b-2。

### 1e. (修正済 2026-10-28) `feedback.f90` の SN イベントログで添字違い (上流由来)
- `feedbk` の `sf_log_properties` ブロック末尾 `typep(ind_part(i))%tag` が、ループ変数 `j` の
  ところで古い `i` を使っていた (上流 `pm/feedback.f90:405` も同じ)。`sf_log_properties=.true.`
  で SN が起きると `typep(0)` を読み、-O2 では segfault、bounds check では
  "Index '0' of dimension 1 of array 'typep'"。`j` に修正。`patch/mizuki/feedback.f90` にも同じ行がある
  (そちらは触っていない)。本番 namelist は `sf_log_properties=.true.` なので本番でも踏み得た。
- 2026-10-29: origin/mizuki acb88172 でログ出力ブロックが MHD のループ注入の前に移動した (こちらにも
  移植済み、hydro ビルドでは bit 一致)。mizuki 側はまだ `typep(ind_part(i))` のまま。小野さんに連絡要。

### 1f. (修正済 2026-11-07) `LIGHT_MPI_COMM=1` で DICE 初期条件が読めない

- to_dice の 4 ファイルに LIGHT 分岐が無い箇所が 7 つあり、`reception(icpu,ilevel)%igrid` / `%up` で
  コンパイルが止まる (`init_flow_fine` ×2, `rho_fine` ×1, `particle_tree` の DICE 追加分 `up`/`maskp` ×4)。
  LIGHT 側は `%pcomm%igrid`、粒子は `reception%pcomm%u(current_property,ind_com(i))` /
  `emission_part(ilevel)%u(current_property, offset_np+ind_com(i)-1)` (添字の順が heavy と逆) で書いた。
- 上流由来: `init_part.f90` の `load_grafic` にある LIGHT 版の粒子再配分 (`emission_part(1)%u` への詰め込み) は
  `sendbuf(icpu)` を増やさず全粒子が同じ場所に上書き、`MPI_ALLTOALL` に渡す送信数が 0、添字が受信側の
  `pcomm%u(1:ncache,1:7)` (列優先) と不一致、`MPI_ISEND` の先頭が `×(twondim+1)` されていない、送るものが無い
  ランクは自分の粒子の詰め直しまで飛ばす、受信側の `reception(icpu,1)%pcomm` ポインタを確保せずに中身を確保
  (level 1 に ghost が無い相手で segfault)、送るものが無いランクで未確保の送信バッファを `deallocate`、の 7 点が
  壊れていた (GitHub `stable` も同じ)。受信側に合わせて書き直した (`to_dice/init_part.f90` と
  `patch/takashi/init_part.f90` 同文)。DICE は `load_dice` (rank 1 が読んで BCAST、各ランクが `cmp_cpumap` で
  選ぶ) なのでこの経路を通らない。
- grafic の検証: `~/ramses_ks_20230630/ics/ics_L1Mpc_64/level_006` (64^3, 1 Mpc, DM+バリオン) を RT=0, DICE 空の
  ビルドで 4 ランク 4 step (levelmin 6, nremap=1)。LIGHT=1 と 0 で `amr/hydro/part/grav` 48 ファイルが bit 一致、
  ランクごとの npart も一致。namelist は scratchpad `run/gh/run.nml` (`&HYDRO_PARAMS` が必須)。
- DICE の検証: `tests/dice_rt_amr.nml` (4 ランク, 6 step, `star_as_nbody=.false.`) で LIGHT=1 と LIGHT=0 の
  `amr/hydro/part/grav/rt` 40 ファイルが bit 一致。LIGHT ビルドは `bin/*.o *.mod` を消してから。

## 2. level ごとの domain decomposition (balance_per_level) — 断念・削除済み (2026-10-29)

細かい level `l` の担当領域と level `l-1` の担当領域が別の場所になる限り、level `l` の oct の親 grid
(とその 3×3×3 近傍) と自分の level `l-1` grid の子 oct がまるごと ghost になり (docs 4c 項 1)、これは
level ごとの均等化そのものと両立しない。本番でも RT subcycle あたりの時間は変わらず、所有者の不整合
(旧 1d) も未解決だった。コードを複雑にしたまま残す価値がないので、分割関係の変更はすべて取り除き
(`load_balance.f90`, `refine_utils.f90`, `flag_formation_sites.f90`, `init_tracer.f90`, `sink_particle.f90` の
patch コピーは削除して上流に戻し、`init_amr/output_amr/read_params/init_part/amr_parameters` は per-level 前の
状態に、`amr_commons` の `bound_key` は 1 次元に、`virtual_boundaries` の `authorize_fine` は元の形に、
`build_comm` の所有者チェックと `debug_missing_owner` は削除)。通信の集約・非同期化 (docs 3 節:
`make_virtual_fine_dp_nv`, `make_virtual_reverse_dp_nv/_start/_end`) は分割と独立なので残した。
出力ファイル形式も元に戻った (per-level 対応コードが書いた `amr_*.out` は末尾に余分なレコードがあるだけ
なので、戻したコードでも restart できる)。設計・検証の記録は `docs/level_domain_decomposition.md`。

## 2a. 単一分割の重みの knob `cost_weight_power` (2026-10-29 追加、本番で効果なし 2026-11-07: 2b 末尾)

待ち (本番の交換 2 本 21% のうち純粋な通信を除く分) は level 15 の grid 数の偏り (max/avg 1.38、
rt - godunov の max/avg 1.31 と一致) から来る。`cost_weighting` はセルの重み `(80+numbp)×niter_cost(l)`
(`niter_cost(l)=Π nsubcycle`、本番では 2^(l-6)) の合計を揃えるだけなので、level の混ざり方が違う
ランク間で level 15 の数は揃わない。`&RUN_PARAMS` の `cost_weight_power` (既定 1.0 = 従来と bit 一致)
で `niter_cost(l)**power` にすると細かい level の重みが支配的になり、level 15 がより均等になる代わりに
grid 総数 (メモリ) の偏りが増える。メモリに余裕があるラン (flat-330638 は最大ランクで 27%) 向け。
`flag1` は integer*4 なので `(80+numbp)×512^power < 2^31` の制限 (power ≲ 1.5、超えると停止メッセージ)。
- 本番での試し方: 同じ restart から `nstepmax` を 3 step 先にして power=1 / 1.3 / 1.5 を流し、
  `Level 15 has N grids (min,max,avg)`、`mem=`、timer 表の `rt - godunov` と交換 2 本の min/max を比べる。
  level 15 の max/avg が 1.1 程度まで下がり、mem の最大が ngridmax に収まれば採用。
- 実装: `patch/takashi/load_balance.f90` (上流コピー + 数行)、`amr_parameters.f90`、`read_params.f90` ×2。

## 2b. RT step のコスト内訳 (2026-10-28 着手)

本番 (896 ランク) の coarse step 時間は level 15 の RT subcycle 1 回あたり 0.39 s で決まり、level 15 を
均等化しても変わらない (docs 4c)。何が 0.39 s なのか (rt godunov の計算か、cooling の hotspot か、
ghost 交換の待ちか) を切り分けるため `rt_step` に細かい timer (`rt - courant / set unew / star feedback /
rev ghostzones / godunov / sources / set uold / ghostzones / upload fine`) を入れた (297d3a22)。
- timer 表は正常終了時 (`clean_end`) にしか出ない。本番で `nstepmax=7578` (restart の nstep_coarse=7575
  から 3 coarse step) にして単一分割で流し、表の min/avg/max を見る。max≫min なら偏り (cooling hotspot
  など)、min も大きければ本当の計算/通信コスト。
- 手元 (8 ランク, levelmax 8, IC の星を光源にして RT 輸送を有効化) では rt godunov が 55%、cooling 11%、
  rev ghostzones 4.6% (min/max 1.5 倍)。1 セル 1 subcycle あたり godunov ≈ 6 µs (bounds check 付き -O1)。
- 手元の `dice_rt_amr.nml` は `star_as_nbody=.true.` で 12 step では星ができないため `rt_advect` が
  立たず RT 輸送が走らない (chemistry だけ)。RT の性能を見るときは `star_as_nbody=.false.` にする。

**本番の timer 表 (flat-330638, 2026-10-29, 単一分割, 896 ランク, 3 coarse step = 1258 s)**:

| 項目 | 平均 s | % | min – max |
|---|---|---|---|
| rt - godunov | 342 | 27 | 260 – 448 |
| rt - ghostzones (forward) | 148 | 12 | 65 – 275 |
| rt - rev ghostzones | 113 | 9 | 21 – 281 |
| cooling | 79 | 6 | 74 – 85 |
| feedback | 177 | 14 | 111 – 288 |
| poisson | 139 | 11 | 138 – 140 |
| hydro - godunov | 131 | 10 | 111 – 150 |
| hydro - rev ghostzones | 62 | 5 | 40 – 95 |

- RT subcycle 数 (3 step): level 15: 1932, 14: 966, 13: 483, 12: 240。grid 数: level 15 2.06M
  (1 ランク 784–3143, 平均 2303), 14 1.71M (729–3494), 13 0.80M (91–3875)。
- 全 level の「セル × subcycle」≈ 4.9e10、godunov 合計 342×896 s → **6.3 µs / セル / subcycle**。
  手元の bounds check 付き -O1 と同じ値なので -O3 で速くなっていない (gather・メモリ律速)。
  `rt_godfine1` 自体の改善 (subcycle ごとの stencil gather と粗 level からの補間のやり直し等) が
  最大の単独項目。
- 通信 2 本 (21%) は純粋な通信 + 遅いランク待ち。godunov の min–max 差 190 s がそのまま待ちに
  乗る。min 側 (rev 21 s, fwd 65 s) から純粋な通信は 1 回 5–17 ms、全体の 7% 前後と推定。
- godunov の偏りは level ごとの grid 数の偏り (`cost_weighting` は全 level の重み付き合計を揃える
  だけ)。cooling は揃っていて hotspot ではない。RT の次は feedback (14%, 2.6 倍の偏り)。
- `output_timer` を拡張し、正常終了時に `timer_ranks.dat` (全 timer のランク別秒数 + level ごとの
  active grid 数) を書くようにした (出力時は `<timing file>.ranks`)。次の本番ジョブで
  godunov の偏りと grid 数の相関、純粋な通信の見積もりを確定する。同じ 3 step を
- per-level で速くならない理由 (2026-10-29): 本番 3 ジョブの coarse step 時間を
  `T = A + B × (level 15 の RT subcycle 数)` で当てはめると単一分割 A=168 s, B=0.39 s、per-level
  A=170 s, B=0.39 s で同一。per-level は level 15 の grid 数を 762–3178 → 2301–2303 に揃え、ghost を
  0.32 → 1.5 (active 比) に増やした。「揃えた利得が ghost 5 倍の通信増で相殺」か「godunov の偏りは
  grid 数ではなく 1 grid あたりの仕事 (補間が要る境界 grid の割合、メモリ配置) の差」かは
  `timer_ranks.dat` の godunov と grid 数の相関で決める。手元 (8 ランク) は元から揃っているので
  判別できない。手元の単一/per-level 比較は IC の星の SN (per-CPU 乱数 `ranf(localseed)`) で
  最初の step から履歴が分かれるので、時間の直接比較には使えない (subcycle 数で割ればセルあたりの
  godunov コストは同じ)。

**rt_godfine1 の削減 (2026-10-29, bit 一致, `patch/takashi/rt_godunov_fine.f90`,
`rt_flux_module.f90`)**: RT は 1 次精度なので 6×6×6 stencil のうち読まれるのは 4×4×4 の
角を除く 32 セルだけ (`cmp_flux_tensors` は元から 32 セルに絞ってあった) なのに、
(1) 3×3×3 の親セル 27 個全部について gather と粗 level からの補間をしていた → 読まれるセルを
含む 7 個 (中心 + 面で接する 6 個) に限定、(2) 面の flux を 27 面/方向計算して 12 面しか
使っていなかった → 12 面に限定、(3) GLF の面 flux を配列一時変数 + 関数呼び出しから成分ごとの
素直なループに (演算順序は同じ)。手元 (8 ランク, dice_rt_amr, star_as_nbody=.false., 12 step) で
`rt - godunov` 16.6 → 5.65 → 5.18 s、全体 39.1 → 25.9 s。hydro/rt/part/grav の出力は全て bit 一致。
本番では rt - godunov が 27% なので全体で 15–20% の短縮見込み。残りは `cmp_flux_tensors`
(32 セル × 7 群の Eddington tensor, 隣接 grid 間で 4 重に計算) と gather。

**rt_godfine1 の内訳と第 2 段 (2026-10-29)**: 上の削減後に `MPI_Wtime` で区間計測 (手元 8 ランク,
rank 1 の累積, 計測用ビルドのみ・コミットしていない): `cmp_flux_tensors` 3.4 s (44%), stencil gather
2.6 s (34%, うち粗 level からの補間 0.3 s), 保存量の更新 1.1 s (14%), GLF 面 flux 0.55 s (7%)。
- `cmp_flux_tensors` を成分ごとに書き下ろし (配列一時変数 `pflux`, `u` と p,q の 2 重ループ、
  セルごとの `if(Np<0) clean_stop` と `cycle` を無くし、`merge` で密度 0・flux 0 を処理)、
  更新ループは 3 方向を 1 回の read-modify-write にまとめた (加える順序は X, Y, Z のまま)。
  手元 (同条件, 単独実行) で `rt - godunov` 5.47 → 4.10 s、全体 27.6 → 24.9 s。
- **bit 一致の確認で注意**: Mac (arm64) の gfortran は既定で a*b+c を FMA に縮約するので、
  数学的に同じ式でも書き方が変わると丸めが変わり、通常ビルドでは出力が bit 一致しない。
  両方を `-ffp-contract=off` でビルドすると Stromgren (4 ランク) と dice_rt_amr (8 ランク, 12 step)
  の出力は bit 一致する (演算順序が同じことの確認はこの方法で行う)。クラスタの ifx も
  `-fp-model source` では FMA 縮約が有効なので、本番でも丸めレベルの差は出る。
- 残りの候補: gather は `rtuold(ncell,nvar)` の配置上セルごとに 28 本のキャッシュラインを
  飛び飛びに読むのが律速で、配列の並びを変えない限り大きくは減らない。tensor の 4 重計算を
  無くすには群ごとにセル当たり 12 個の tensor を全セル分持つ必要があり (メモリ +40%)、
  gather が増えるので得にならない見込み。次は feedback (本番 14%, 2.6 倍のばらつき)。

**feedback timer の正体と star_formation の削減 (2026-10-29, bit 一致, `patch/takashi/star_formation.f90`,
`amr_step.f90`)**: `feedback` timer は `thermal_feedback` と `star_formation` の両方を含んでいた。
`star_formation` を `star formation` timer に分けると手元 (8 ランク, dice_rt_amr, 12 step) では
feedback 0.10 s / star formation 1.35 s なので、本番の 14% (2.6 倍のばらつき) はほぼ star_formation。
原因は sf_model=6 の経路で葉セル 1 個ごとに `getnbor` と `get_27_neighbors` (それぞれ 1 セルの
`get3cubefather`) を密度判定の前に呼んでいたこと。(1) sf_model 6/7 かつ sf_virial=.false. では
密度判定 `d<=d0` を先に適用、(2) これらのモデルで使わない `getnbor` (面近傍 d1–d6) を省略、
(3) 27 近傍を grid ベクトル単位で求める `get_27_neighbors_grid` (grid の 8 セルは同じ 3×3×3 親セルを
共有するので `get3cubefather` は nvector grid に 1 回、近傍の対応は `get_27_neighbors` と同一) を追加。
手元で star formation 1.35 → 0.22 s、出力は bit 一致 (近傍の添字と和の順序が同じなので FMA の影響なし)。
sf_virial=.true. や sf_model<6 の経路は元のまま。冒頭の `uold(1:4)` の `make_virtual_fine_dp` 4 回も
`make_virtual_fine_dp_nv(uold,nvar,1,4,ilevel)` 1 回にまとめた (bit 一致、手元 0.22 → 0.21 s)。
本番ではこの交換が level ごと・fine step ごとのブロッキングなので、ばらつきの待ちはここに残る。

**RT の逆方向交換を光子密度 7 変数に限定 (2026-10-29, bit 一致, `amr_step.f90`, `virtual_boundaries.f90`,
`amr_commons.f90`)**: RT subcycle の `make_virtual_reverse_dp_start/end(rtunew)` は 28 変数全部を
送っていたが、仮想セルの rtunew は `rt_set_unew` で 0 にされ、書き込むのは `star_RT_vsweep`
(`rtunew(indp,iGroups(ip))`、光子密度のみ) だけで、`rt_godunov_fine`/`add_rt_sources` は仮想セルを
触らない。そこで start/end に変数の刻み `ivstep` を追加し (旧呼び出しは 1)、RT は
`(rtunew,nrtvar,1,iGroups(nGroups),ndim+1,...)` で密度だけ交換する (データ 1/4)。手元で
`rt - rev ghostzones` 1.32 → 1.08 s、出力は bit 一致 (flux 変数への 0 加算が無くなるだけ)。
本番では純粋な通信分 (min 側 21 s/3 step) が減る。待ちの分は変わらない。
注意: `amr_commons` の型を変えたので、`bin/` の `.o`/`.mod` を全部消してフルビルドが必要
(Makefile にモジュール依存が無く、古い `.mod` と不整合になる)。

**本番の再計測 (flat-334557 / 334558, 2026-11-07, 896 ランク, 同じ restart, 3 coarse step, ログは
`patch/takashi/logs/`)**: f98dad81 のバイナリ (RT kernel 削減 + reverse 7 変数 + forward/uold の重ね合わせ +
star_formation の gather 削減) で、334557 は `cost_weight_power=1.0`、334558 は 1.3。

| 項目 | flat-330638 (10/29) | 334557 (power 1.0) | 334558 (power 1.3) |
|---|---|---|---|
| 3 step 合計 | 1258 s | 692 s | 702 s |
| level 15 RT subcycle 数 | 1932 | 1199 | 1257 |
| rt - godunov 平均 (min–max) | 342 (260–448) | 52 (46–64) | 55 (49–67) |
| rt - ghostzones | 148 (65–275) | 45 (25–68) | 48 (26–76) |
| rt - rev ghostzones | 113 (21–281) | 22 (6–42) | 24 (6–41) |
| cooling | 79 (74–85) | 54 (51–65) | 57 (52–68) |
| star formation (+feedback) | 177 (111–288) | 90+3 (58–133) | 89+3 (59–127) |
| poisson | 139 | 148 | 148 |
| hydro - godunov | 131 (111–150) | 135 (117–155) | 135 (120–155) |
| level 15 grids min/max (最終) | 784 / 3143 | 766 / 3259 | 1244 / 3129 |

- 物理が分岐するので subcycle 数が違う。`T = A + B × subcycle` の当てはめ: 旧 A≈155 s, B≈0.41 s/回 →
  新 A≈140–170 s, B≈0.19–0.21 s/回。**RT subcycle 1 回のコストが半分**、切片はほぼ同じ。同じ subcycle 数
  なら旧 957 s 相当が 692 s (28% 短縮)。
- rt - godunov は subcycle あたり 0.177 → 0.044 s (4 倍)。rt - ghostzones の min は 65 → 25 s、rev の min は
  21 → 6 s (subcycle あたりでもそれぞれ 2 倍, 1.6 倍)。star formation は min 111 → 58 s (fine step 数は同じ)。
- 残りの順位 (power 1.0): poisson 21%、hydro - godunov 19.5%、star formation 13%、hydro - rev ghostzones 10%、
  cooling 8%、rt - godunov 7.6%、rt - ghostzones 6.5%。RT は合計 18% まで下がり、次は poisson と hydro。
- `cost_weight_power=1.3` の効果: level 15 の min が 766 → 1244 に上がるが max は 3259 → 3129 (4%) しか
  下がらず、step 2 (両者とも 512 subcycle) の時間は 241.2 s と 242.1 s で同じ。rt - godunov の max も
  64 → 67 で変わらない。**採用しない** (既定 1.0 のまま)。step 1 は restart ファイルの分割なので両者同一。

**ランク別 timer の解析 (`patch/takashi/logs/timer_ranks_334558.dat`, power 1.3 のラン, 2026-11-07)**:
- rt - godunov (49–67 s)、cooling (52–68 s)、hydro - godunov (120–155 s) は**同じランクで同じ比率で遅い**
  (相関 0.82 / 0.73 / 0.54)。grid 数との相関は弱く (level 15 とは **負**の相関 -0.45、subcycle 数で重み付けた
  仕事量とも 0.1–0.5)、level 15 を最も多く持つ rank 817 (3040 grids) が最速、level 13/14 を多く持つ
  rank 770/829 (L15 1700–1850, L14 2900–3400) が最遅。grid-subcycle あたりの rt - godunov は 10.5–14.2 µs
  (35% の幅)。ノード単位 (112/128 ランクのブロック) の分散は 3% なのでハードの差ではない。
  → 偏りは「grid 数」ではなく「grid あたりの仕事」。cooling (セル局所) まで同じ比率で動くので、
  refinement 境界の補間だけでは説明できず、ランク内の grid のメモリ配置 (load balance で受け取った grid が
  free list に散らばる) によるキャッシュ効率の差が最有力。確かめるなら load balance 直後に grid を
  Hilbert 順に並べ直す (全配列の付け替え、大仕事) か、`timer_ranks` の「最後の load balance で受け取った
  grid 数」(`recv level l` / `recv total`、2026-11-07 追加、`nrecv_lb` in `amr_commons`) と grid あたりの
  godunov 時間の相関を見る。
- star formation (59–127 s) は grid 数と無相関 (0.1)、rho と 0.55、refine と -0.56。粒子数の行が無かったので
  判定できず、`timer_ranks` に `particles` / `stars` / `ghost level l` の行を足した (次のランから出る)。
- `timer_ranks.dat` の名前は `timer_ranks_<PBS_JOBID の数字>.dat` (SLURM なら `SLURM_JOB_ID`、どちらも無ければ
  `timer_ranks_n<nstep_coarse>.dat`) に変えた。手元 (`PBS_JOBID=123456.flat`) で確認済み。

**追加行つきの再計測 (flat-337110 / `timer_ranks_337110.dat`, power 1.0, 2026-11-07, 677 s; 334557 と同じ設定)**:
- `recv total` は復帰後の全 load balance の累積 (復帰直後の 1 回目で約 360 万 grid が動く)、
  `recv level l` は最後の 1 回だけ (全 grid の 4.3%)。ローカルの DICE テストでも同じ比率なので仕様どおり。
- **「受け取った grid のメモリ配置」説は否定**: grid-subcycle あたりの rt - godunov は、最後の load balance で
  受け取った割合の四分位ごとに 12.23 / 12.05 / 12.06 / 12.09 µs で差がない。累積 `recv total` とも相関 -0.15。
  ghost 割合とも 0.04、粒子数とも -0.06。星の数とは 0.24 (四分位で 12.0→12.3 µs、3%) で主因ではない。
- **grid 数ではまったく説明できない**: level 11–15 の grid 数での最小二乗は R² ≤ 0.1 (rt - godunov、cooling、
  hydro - godunov とも)。cost weighting は仕事量 (4.60–4.64M) をよく揃えているのに、計算時間は
  rt - godunov 50.8–69.2 s、cooling 53.7–68.8 s、hydro - godunov 91–123 s と最大 +20% ばらつく。
  3 つのタイマーは純粋な計算部分 (非同期通信の end 待ちは `rt - rev ghostzones` / `rt - ghostzones` /
  `hydro - rev ghostzones` に入る) なので通信待ちの混入ではない。
- **遅いランクは別ジョブでも同じ**: 334558 (power 1.3、分割が違う) と 337110 で grid あたりコストの相関は
  rt - godunov 0.72、cooling 0.82、star formation 0.83。平均より 10% 以上遅いランク (41, 80, 81, 83, 109,
  407, 645 …) は両方に出る。ランク番号の mod 2/4/8/16 (コア位置) には依存せず、112/128 ブロック (ノード) で
  説明できる分散は 10–15%。ログにノード名が無いので「同じノードに割り当てられた」可能性は否定できないが
  (2026-11-07 から起動時に `Ranks i - j on <node>` を出し、`timer_ranks_*.dat` の 2 行目 `# host` に
  ランクごとのノード名が入るので、次回以降のラン同士なら判別できる。334558/337110 については
  `qstat -xf <jobid>` の `exec_host` が残っていればそれで)、
  Hilbert 順でランク番号 ≈ 空間領域なので、領域 (中心部の高密度ガス?) 依存の grid あたりコストと見るのが自然。
  rank 番号での自己相関 (lag 1–2 で 0.3–0.4) も領域説と整合する。
- 同期で吸収される待ち時間: star formation 66–145 s、rt - ghostzones 21–72 s、hydro - rev ghostzones 45–88 s。
  合計の最小 132 s 対平均 213 s、速いランクは遅いランクより ~80 s 余計に待つ (677 s の 12%)。
- **次の手**: grid 数ではなく実測時間で重みを補正する。各ランクの直前 coarse step の計算時間 (rt - godunov +
  cooling + hydro - godunov) / 仕事量 を係数にして、そのランクが持つ grid の cost weight に掛けてから
  `cmp_new_cpu_map` に渡す (ランク ≈ 領域なので領域依存コストもハード差も両方拾える)。期待値は最大 ~10%。
  領域依存の原因そのもの (cooling の反復回数? Riemann solver? 高密度での何か) は `rt_solve_cooling` の反復回数を
  ランク別に数えれば分かる。

## 2c. 実測時間で cost weight を補正する `cost_time_weighting` (2026-11-07 追加、本番未検証)

2b の結論 (grid 数では説明できない、ランクごとに再現する +20% の計算時間差) への対処。
`&RUN_PARAMS` に `cost_time_weighting=.true.` (既定 `.false.`) と `cost_time_relax` (既定 1.0)。

- `cmp_new_cpu_map` (`load_balance.f90`) で、各ランクが自分の leaf cell の cost (`flag1`、level 重み × (80+粒子数))
  を集計した直後に、`cost_timer_total` (`update_time.f90`、`hydro - godunov` + `cooling` + `rt - godunov` の
  timer 合計) の前回 load balance からの増分 / 自分の cost 合計 を「単位 cost あたりの実測時間」とし、
  測定できたランクの平均で割って係数 `cost_rank_fac` (`amr_commons`) にする。[0.5, 2] に clip、
  `cost_time_relax` で前回の係数と混ぜる (1.0 なら最新値のみ)。
- 係数は `cost_loc` (ランクの総 cost) と境界探索の `cell_cost` の両方に掛ける (`flag1` は整数のまま)。
  同じランクの cell は全部同じ係数なので「そのランクの領域は単位 cost あたり X 倍遅い」という補正。
  領域依存でもハード依存でも、ランク ≈ 領域 (Hilbert 順) で動く grid は 4% 程度なので両方拾える。
- 最初の load balance (初期条件、restart 直後) は計測区間が無いので係数 1 (`cost_time_valid`)。
  `cost_time_weighting=.false.` では `cost_rank_fac=1` を掛けるだけで出力は bit 一致
  (手元 `tests/dice_rt_amr.nml` 相当、4 ランク 6 step: 変更前バイナリと 40 ファイル cmp 一致。
  `.true.` では係数 0.98–1.02 で分割がわずかに変わり、`hydro_stats.py` の総和は丸め誤差の範囲で一致、
  セル数は同じ)。
- 出力: rank 1 が毎回 `Load balance: time cost factor min= max= (N ranks measured)` を出す。
  `timer_ranks_*.dat` に `cost factor` 行 (各ランクの係数)。
- 本番で確認すること: 係数の min/max が 0.8–1.2 程度に収まるか、振動しないか (収まらなければ
  `cost_time_relax=0.5`)、計算 3 タイマーの max/avg が 1.2 → 1.05 程度に縮むか、
  待ち (star formation / rt - ghostzones / hydro - rev ghostzones) の合計が減るか。
  memory (grid 数) は不均等になる方向なので `mem=` の最大値も見る。

**本番での結果 (flat-337199 / `timer_ranks_337199.dat`, `cost_time_weighting=.true.`, 2026-11-08, 725 s)**:
- 係数は 0.93–1.08 (step 2)、0.93–1.07 (step 3)。分割の cost は 13.7% の幅になったが、計算 3 タイマーの
  max/avg は 1.21 → 1.18 とほとんど縮まず、遅いランクは前と同じ (41, 80, 81, 83, 109, 352, 407, 645, 829, 871)。
  rank 109 (累積で最遅、1.16–1.21 倍) の係数は 0.950。**係数は遅いランクを捉えていない**。
- 手元 (4 ランク) では係数と累積の sum3/work が 0.996 vs 0.994 のように一致するので測定機構は正しい。
  → 本番の遅さは coarse step ごとに**間欠的** (rank 109 は step 2 では遅くなく、1 と 3 で ~1.3 倍) で、
  しかもジョブをまたいで同じランク・同じ step で再現する (同じ snapshot からの restart なので当然)。
- 遅いランクは**単独** (隣の 108/110 は正常) で、cooling・rt - godunov・hydro - godunov が**揃って** ~15% 遅い。
  ノードは今回 cna037–cna201 で、遅いランクのノード内位置はばらばら、ノード平均の差は ±4% → ハードではない。
  手元の出力には非正規化数は無い (RT は 0 か正規数) ので denormal 説も根拠薄。
- 残る候補: そのランクの領域で何かが step によって起きる (SN 直後の高温セルで cooling の反復が増える、
  refinement 境界の補間が増える等) が、3 カーネル一律というのが説明しにくい。まず **step ごと・ランクごと**の
  時間を取る: `&RUN_PARAMS timer_ranks_every=1` で coarse step ごとに `timer_ranks_<jobid>_n<step>.dat`
  (累積値) を書くようにした (2026-11-08)。差分で step ごとの時間が出る。level ごとの内訳が要るなら
  godunov/cooling の timer を level 別に分ける。SN 説 (岡本さん: SN が起きた所で遅くなるのは十分ありえる) を
  直接見るため、ランクが注入した SN 数の累積 `nsn_rank` (`amr_commons`、`feedback.f90` で加算) を
  `SN events` 行として同じファイルに出す。step ごとの差分と cooling / godunov の差分を突き合わせる。
  `sf_log_properties` の `stars_*.out` は時刻が無いので step に分けられない。
  なお `stars` 行は粒子スロット (全スロット − free list) から数えるように直した: coarse step 終了時点では
  `merge_tree_fine` 後で子 grid の粒子リストが古いままなので、grid リスト経由だと 3 倍に数えていた
  (終了時・出力時のファイルは正しかった)。手元で step ごとの `stars` が終了時と一致することを確認。
- 本番のノード: 28 ランク/ノード × 32 ノード。

**解決 (flat-337718 / `timer_ranks_337718_n7576..7578.dat`, `timer_ranks_every=1`, 2026-11-09, 672 s)**:
- 遅いランク (41, 80, 81, 83, 109, 352, 407, 645, 829, 871 …) が遅いのは **step 1 だけ** (hydro - godunov で
  1.3–1.6 倍)。step 2・3 では 0.96–1.07 で平均並み。SN は 3 step 合計 3 個で無関係。
- 原因は **restart 直後の分割**: 最初の coarse step では load balance がスキップされ (amr_step の
  "Skip first load balance because it has been performed before file dump")、snapshot の cpu map のまま走る。
  その分割の仕事量 (grid × 2^(level−7)) は max/avg **1.68**、L15 grid 数 789–4058。rank 109 は仕事量 1.68 倍で
  時間 1.54 倍、仕事量あたりでは 0.92 と普通。step 2 の最初の load balance で 324 万 grid が動いて
  max/avg 1.09、step 3 で 1.03 になる。snapshot を作ったランは今と cost 模型が違ったのだろう。
  2b の「grid あたりコストがランクで違う」「間欠的」は全部これの見間違い (累積値しか無かったため)。
- 定常状態 (step 2・3) の計算 3 タイマーの仕事量あたりの偏りは max/avg 1.13–1.14、遅いのは連続した
  683–689 (1.03–1.14) で領域依存。ここは `cost_time_weighting` が効きうる範囲 (最大 ~10%)。
- 337199 (`cost_time_weighting=.true.`) の step 3 が 52 s 遅かったのは分割のせいではなく、物理が僅かに
  ずれて dt が変わり L15 の RT subcycle が 280 → 512 になったため (337110 と 337718 は dt まで同一)。
- **3 step の restart テストの正しいやり方**: `&RUN_PARAMS restart_remap=.true.` で最初の step から
  load balance させる (step 1 の 249 s のうち ~60 s は分割の偏り)。または step 1 を比較から外す。
  `cost_time_weighting` を試すなら、係数が restart 分割で測られないように `restart_remap=.true.` と組で。

## 2d. OpenMP 試作: `godunov_fine` だけスレッド化 (2026-11-09)

背景: 本番は 112 コア × 32 ノードで 896 ランク (28 ランク/ノード) しか使っていない (3584 ランクでは
10 倍以上遅い)。空いている 3/4 のコアをスレッドで使うための第一歩として、まず `godunov_fine` だけを
OpenMP 化し、1 ノードで 1・2・3 スレッドの伸びとメモリ帯域の余裕を測る。RAMSES-yOMP
(Han+2025, github.com/sanhancluster/RAMSES-yOMP) を参考にしたが、あちらは rt/ が手つかずで
上流とも大きく離れているので、指示文の入れ方だけ真似て patch 側に自前で入れる。

- Makefile に `OPENMP = 0/1` (GNU `-fopenmp`, Intel `-qopenmp`; `bin/Makefile` と `Makefile.xd_intel`)。
- `godunov_fine` の vector sweep ループを `!$omp parallel do` (dynamic,2)。`godfine1` は自分の grid の
  セルにしか書かないが、粗い隣接セルへの flux 補正 (`unew/divu/enew(ind_buffer)`) は複数の fine grid が
  同じ粗セルに足すので `!$omp atomic update` (6 箇所)。足す順が変わるのでスレッド数 >1 では丸め誤差
  レベルで結果が変わる (1 スレッドは MPI 版と bit 一致)。
- `save` 作業配列は thread-unsafe なので除去: `godunov_fine.f90` (23)、呼ばれる `umuscl.f90` (9)、
  `godunov_utils.f90` (24)、`interpol_hydro.f90` (11)、`nbors_utils.f90` (10) を patch/takashi にコピーして
  `,save` を外した。`godfine1` の `gloc/ploc/req_loc/peq_loc/req2/peq2` は save+初期化で 0 を保証していた
  ので毎回明示的に 0 にする。`nbors_utils` の `iii/jjj/kkk` は定数表 (読むだけ) なので初期化付きのまま。
- スタック: `save` を外すと `-fopenmp` (gfortran は `-frecursive` 相当) で作業配列がスタックに乗る。
  godfine1 + unsplit で 1 スレッドあたり ~8 MB 必要。**`OMP_STACKSIZE=128m` と `ulimit -s`
  (main スレッド) を必ず設定**。Intel なら `-heap-arrays` も選択肢。
- OPENMP=0 のビルドは `save` 除去前と bit 一致であること (gfortran は大きい局所配列を自動的に静的領域に
  置くので挙動は同じ)。

**結果 (手元 Mac, M 系 6P+6E コア, gfortran -O3, DICE テスト levelmax 8, 6 coarse step, 2026-11-09)**:
- OPENMP=0 ビルド: 参照 (save 除去前) と 40 ファイル bit 一致。OPENMP=1 の 1 スレッドも bit 一致。
  2・3 スレッドは hydro/rt の出力が丸め誤差で変わる (`hydro_stats.py` の総和・min/max の相対差 ≤ 1.3e-10)。
- 1 ランク × スレッド数 (P コア内) の `hydro - godunov`: 1: 10.04 s、2: 5.53 s (1.81 倍)、4: 2.98 s (**3.37 倍**)。
  godunov は 4 スレッドまでほぼ素直に伸び、メモリ帯域は律速していない。
- 4 ランク × 3 スレッド (12 スレッド、E コア込み) では 3.14 → 1.98 s (1.59 倍) にとどまる: E コアのせいで、
  スレッド化の問題ではない。
- 次: 同じやり方で `rt_godunov_fine` (`rt_godunov_utils`/`rt_interpol_hydro` の save 除去込み) と
  `cooling_fine_part` (セル局所、`rt_solve_cooling` は save を使う? 要確認)、multigrid の smoothing、
  `virtual_boundaries` の pack/unpack、`rho_fine` (スレッド別バッファ)。クラスタでは Intel `-qopenmp` で
  28 ランク × 4 スレッド、`OMP_STACKSIZE=128m`、`ulimit -s unlimited`、`I_MPI_PIN_DOMAIN=omp` を job script に。
  cooling の `rt_solve_cooling` は乱数を使わないので決定的にできるはず。
- 未追跡ファイル: `patch/takashi/{umuscl,godunov_utils,interpol_hydro,nbors_utils}.f90` (コミット時に git add)。
- `OMP_STACKSIZE` を指定しないとどうなるか (手元、1 ランク × 2 スレッド): 最初の godunov でワーカースレッドの
  スタックが溢れる。macOS では stack overflow が SIGBUS (= signal 10) になり、RAMSES が SIGUSR1 (Linux の 10)
  用に `SIGNAL(10,output_signal)` で掛けているハンドラがそれを捕まえて "SIGNAL 10: Output will be written..."
  を出して戻る → 同じ命令で再び溢れる → **無限ループでログが 14 GB** になった (kill して削除)。Linux では
  SIGSEGV (11) なので普通に落ちる。いずれにせよ `OMP_STACKSIZE` は必須。
- `jobfiles/run_hybrid_intel.rb`: `run_partial_intel.rb` の hybrid 版 (`--cpus-per-task=4`, `SRUN_CPUS_PER_TASK`,
  `OMP_NUM_THREADS=4`, `OMP_STACKSIZE=64m`, `OMP_PLACES=cores OMP_PROC_BIND=close`, `srun --cpu-bind=cores`,
  初回確認用に `I_MPI_DEBUG=4` と `--cpu-bind=verbose` (`PinDebug`))。

### 2d-2. rt_godunov_fine と cooling_fine_part の OpenMP 化 (2026-11-09)

- `rt_godunov_fine`: sweep ループを `parallel do`、粗セルへの `rtunew(ind_buffer)` 補正 2 箇所を atomic。
  `rt_godfine1` と `rt_interpol_hydro.f90` (rt/ から patch にコピー)、`rt_flux_module.f90` の
  `cmp_rt_faces` の `save` 作業配列を除去。`rt_per_bnd` (save + .false. 初期化) は毎回 .false. に。
- `cooling_fine` / `cooling_fine_part`: sweep ループを `parallel do` (part の `do while` は `do igrid=igrid0,nlist,nvector` に)。
  `coolfine1` の `save` 配列を除去。`cooling_on`、`Np_boost`、`dNpdt`、`Fp_boost`、`dFpdt`、`polytropic_constant`、
  `err`、`emag` は「条件付きでしか代入されないのに無条件に読まれる」ので毎回明示的に初期化
  (元のコードは save の静的領域がゼロ初期化されることに暗黙に頼っていた)。`is_emitted`/`list` の
  save は意図的な永続データなのでそのまま (並列ループの外でしか書かない)。
- `rt_cooling_module`: `cool_step` の `save` スカラー約 50 個を除去し、物理が off のとき代入されないもの
  (`metal_prime` は一度も代入されずに読まれていた、`dustAbs` 等) を毎回 0 に。`getMu` の save も除去。
  統計カウンタ `n_cool_cells`/`tot_cool_loopcnt`/`loopCodes`/`max_cool_loopcnt` は `!$omp atomic update`
  (`max` も atomic で可)。初回表示の `ifirst` ブロックは `critical`。
- **bit 一致の確認で FMA に嵌った**: OPENMP=0 ビルドが変更前と 20 ファイル違った (相対 1e-15〜1e-10)。
  `save` を外すと配列が静的領域からスタックに移ってコード生成が変わり、gfortran 既定の
  `-ffp-contract=fast` で FMA にまとめる式が変わるため。旧・新コードを `-ffp-contract=off` で
  ビルドすると 40 ファイル bit 一致 → 論理は同じ。Makefile に `BITREPRO = 0/1` (GNU `-ffp-contract=off`、
  Intel `-fp-model strict`) を足した。**以後の bit 比較は BITREPRO=1 のビルド同士で行う**
  (本番ビルドでは不要)。Apple Silicon も Xeon も FMA 命令を持つので両方で起きる。

- **競合の発見**: cooling だけ 3 スレッドにすると結果が 1e-1 ずれた (RT だけなら bit 一致、hydro だけなら 7e-14)。
  原因は `coolrates_module.f90` の `inp_coolrates_table` が `save` ローカルで「前回の T と補間の添字・重み」を
  キャッシュしていたこと (`Tlast`, `iT`, `yy`…)。スレッド間で共有されると別スレッドの T の係数を使う。
  `!$omp threadprivate(...)` でスレッドごとのキャッシュにして解決 (キャッシュの効果はそのまま)。
  ほかに `compCoolrate` の save スカラー、`cmp_chem_eq`/`rt_cmp_metals` の初期化付きローカル (暗黙の save) も除去。
  **教訓: `save` だけでなく宣言時初期化 (`real(dp)::x=0`) も暗黙の save でスレッド共有になる。**
  チェック方法: `grep -n "::.*=" | grep -v parameter` と `grep -in save`。
  修正後: cooling のみ 3 スレッドでシリアルと 40 ファイル bit 一致、全成分 3 スレッドで総和の相対差 7e-14
  (hydro の粗セル補正の atomic 順序)、1 スレッドは bit 一致。
- `&RUN_PARAMS omp_hydro / omp_rt / omp_cooling` (既定 .true.) で成分ごとに `parallel do` を `if(...)` で
  切れる。この切り分けテスト (1 成分ずつ 3 スレッド、シリアルと cmp) が競合の発見に効いた。
- **スレッドの伸び** (1 ランク、P コア、`-ffp-contract=off`、DICE levelmax 8、6 step):

  | スレッド | hydro - godunov | rt - godunov | cooling | TOTAL |
  |---|---|---|---|---|
  | 1 | 9.51 s | 15.20 s | 25.34 s | 64.7 s |
  | 2 | 5.28 (1.80×) | 7.99 (1.90×) | 13.99 (1.81×) | 41.8 s |
  | 4 | 3.05 (3.12×) | 4.36 (3.48×) | 8.11 (3.12×) | 30.1 s (2.15×) |

  3 カーネルとも 4 スレッドで 3.1–3.5 倍。TOTAL の残りはスレッド化していない部分 (poisson, 粒子, refine…)。

- **threadprivate と macOS**: 上の cooling の数字は競合があった版のもの。`inp_coolrates_table` のキャッシュを
  `threadprivate` にして直すと 4 スレッドで cooling 12.85 s (1.97 倍) に落ちた。macOS の gfortran は TLS を
  エミュレーション (`__emutls_get_address`、変数アクセスごとに関数呼び出し) で実装しているため。
  キャッシュを外して毎回 log10 を計算する版は 1 スレッド 29.6 s (キャッシュ有り 25.3 s、+17%)、
  4 スレッド 8.27 s (3.57 倍)。Linux/Intel の TLS は数命令なので本番は threadprivate キャッシュが既定、
  Makefile `COOLRATES_NOCACHE = 1` (`-DCOOLRATES_NOCACHE`) でキャッシュ無しを選べる (Mac でのスレッド計測用)。
  **Mac で cooling のスレッド伸びを測るときは `COOLRATES_NOCACHE=1` で。**
- cooling の `schedule(dynamic,1)` (32 grid 単位) は `dynamic,2` より 4 スレッドで 3.5% 速い (8.27 → 7.98 s、
  NOCACHE 版)。セルごとの反復回数の差が大きいので細かい単位にした。hydro/RT は `dynamic,2` のまま。

### 2d-3. hybrid の本番テスト 1 回目 (hybrd-339635、2026-11-11)

896 ランク × `OMP_NUM_THREADS=4`、32 ノード、`run_hybrid_intel.rb` (`--cpus-per-task=4`)、`cg_auto=.true.`、
同じ snapshot からの 3 step。比較対象はスレッド無しの 339614 (同じ namelist)。

- 全体: 674 → 473 s (1.42 倍)。coarse step ごとに 283/216/166 → 198/154/120 s。RT subcycle 数と dt の履歴は
  2 ランで同じ (L15 subcycle 656/510/256、dt<1e-5 の fine step 0/2/246 対 0/2/256) なので直接比べられる。
- **ただし速くなったのはスレッド化した部分だけではない**: 1 step あたりの平均で
  poisson 30.5 → 18.2 s (1.68 倍、スレッド化していない)、hydro-rev ghostzones 1.25 倍、star formation 1.2–1.4 倍、
  refine 1.3 倍。スレッド化した hydro-godunov は 35 → 25 s (1.38 倍)、rt-godunov 1.51 倍、cooling 1.53 倍で、
  スレッド無しの poisson より伸びが小さい。
- 原因は CPU の割り当て。`--cpu-bind=verbose` のマスクは task 0 が論理 CPU {0,1,112,113}、task 1 が {56,57,168,169}
  で、**ノードは 224 論理 CPU (112 コア × HT 2)、Slurm の `--cpus-per-task=4` は論理 CPU 4 個 = 物理コア 2 個**
  だった。4 スレッドは 2 コアを HT で共有していた (`OMP_PLACES=cores` なので 1 place に 2 スレッド)。
  `sacct` では両ジョブともノード専有 (32 × 224 = 7168 CPU) で、flat 339614 の srun step は 896 task に 7168 CPU
  (1 task 8 論理 CPU = 4 コア)、hybrid 339635 の step は 3584 CPU (1 task 4 論理 CPU = 2 コア)。つまり **flat の方が
  1 ランクあたりのコアは多かった** のに、スレッド化していない poisson まで hybrid の方が 1.68 倍速い。
  flat script (`run_partial_intel.rb`) には `--cpu-bind` が無いので、どう束縛されたか (Slurm 既定か Intel MPI の
  pinning か、HT の兄弟に 2 ランクが乗っていなかったか) はログから分からない。hydro-godunov のような局所計算まで
  1.38 倍速いので、flat ではコアの共有か周波数の問題があった疑いが強い。
  いずれにせよ 339614 と 339635 の差は OpenMP の効果ではなく、hybrid 内の 4 スレッドは 2 コアの HT 共有なので
  スレッドの正味の効果はこのランからは出せない。
- **flat の束縛を見る方法**: `run_partial_intel.rb` の srun に `--cpu-bind=verbose` を足す (型を指定しなければ
  既定の方針のままマスクだけ出る)。
- **次のテスト**: `run_hybrid_intel.rb` を `ThreadsPerCore=2`, `CoresPerTask=4` (→ `--cpus-per-task=8`) に直した
  (2026-11-11)。これで 1 task に物理コア 4 個、`OMP_PLACES=cores` で 1 コア 1 スレッド。
  (a) `NumThreadsPerProc=4` の hybrid、(b) 同じ script で `NumThreadsPerProc=1` の flat (同じ配置での基準) の 2 本を流し、
  (b) に対する (a) の hydro-godunov/rt-godunov/cooling の比が OpenMP の正味の効果 (手元では 3.1–3.7 倍)。
  (b) と 339614 の差が「配置だけの効果」。

### 2e. poisson まわりの上流変更の扱い (2026-11-09 時点、保留)

`poisson/` は 2023-09 の上流のまま。上流 (ramses-organisation) のその後の poisson 変更は 3 つで、PCG の類は無い。
- #313 (2026-09) multigrid の `safe_mode` を 10 回の解で自動解除し restart に保存。我々のコードでは一度入ると
  そのランの間固着し (coarse MG の smoothing が 2 倍 → 交換も 2 倍)、入ったことはログに出ない。ただし
  restart で必ずリセットされるので固着は最長 1 セグメント (24 h)。3 step テストでは全 level `Step= 4` で
  収束していて trigger の気配は無い。**24 時間ランを再開するときに**、`multigrid_fine_commons.f90` を patch に
  コピーして safe mode の on/off を rank 1 が 1 行出す診断を入れ、1 セグメント見てから countdown を移植するか決める。
- #345 (2026-10) multigrid 通信で `flag2` を作業領域に使い回す稀なバグ (別 level の grid index を有効化、
  または `-mod(flag2,ngridmax)` が 0 になって `lookup_mg(0)` で落ちる; issue #344、孤立銀河で 23958 step 目に発生)
  の修正。**2026-11-09 に取り込み済み**: `poisson_commons.f90`, `multigrid_fine_commons.f90`,
  `multigrid_fine_fine.f90` を patch にコピーして `lookup_mg_backup(1:ngridmax)` に退避 (`init_amr.f90` で確保)、
  scan flag は `flag2` にそのまま 0/1 を入れる。通常ランでは結果は変わらない (BITREPRO=1 で変更前と bit 一致を確認)。
- #137 (2025-06) CG のバグ修正は 2025-04 の最適化で入ったバグなので我々には無関係。`cg_levelmin` テストはそのまま可。

### 2f. `cg_levelmin` テスト (flat-338146: cg_levelmin=13、flat-338148: cg_levelmin=11、2026-11-11)

同じバイナリ (9321a939)、896 ランク、337718 と同じ snapshot から 3 coarse step。`restart_remap` は無しなので
step 1 は snapshot の分割 (work max/avg 1.69) のまま (2b と同じ)。

- **poisson の時間 (1 step あたり、全ランクほぼ同じ)**: MG のみ (337718) 55.9 / 48.4 / 47.3 s →
  cg_levelmin=13 で 34.6 / 30.6 / 30.6 s → cg_levelmin=11 で 34.7 / 29.9 / 29.6 s。
  **CG を level 13–15 に使うと poisson が 36% (約 17 s/step、全体の 9%) 減る**。11 まで下げると
  さらに 0.7–1.0 s/step (2–3%) 短いが、ラン内のランク間・step 間のばらつき (約 1 s) と同程度で 1 ランずつでは
  判別できない (level 11/12 は grid 数が少なく、CG 105/83 反復と MG の V-cycle 4 回で手間がほぼ同じ)。
  **`cg_levelmin=13` を採用** (11 でも害はない; 差を確かめるなら複数ラン必要)。
  反復数: L13 73–74、L14 67–68、L15 52 (epsilon は同じ 1e-4)。
- TOTAL (3 step) は 674 → 672 → 687 s で変わっていないが、これは poisson の得が消えたのではなく、
  3 ランの時間発展が最初の SN の後で分岐して RT subcycle 数と cooling の反復が変わったため:
  L15 の RT subcycle 合計 1199 (337718) / 1452 (338146) / 1518 (338148)、rt-godunov 55.8 / 67.6 / 70.8 s、
  cooling 57.7 / 69.8 / 73.0 s。hydro-godunov は 103.6 / 103.9 / 103.7 s で同じ (fine step 数 769 も同じ)。
  **ランごとの比較は poisson の行で見ること**。ノード集合も 3 ランで違う (hydro-godunov が同じなのでノード速度の差は無い)。
- 副産物: 337718 では fine step 12 で SN が 0.17 Msun のセルに入って T2 = 2.3e10 K になり、dt が
  5.4e-5 → 3.0e-6 (18 分の 1) に落ちて、その後 1 step 1.1 倍ずつしか戻らない (769 fine step のうち 265 が dt < 1e-5)。
  この coarse step 群が進めた時間は 0.0131 (337718) に対し 0.0175 (338146): 同じ計算量で 34% 違う。
  **SN が軽いセルに入ったときの dt 崩壊が実効速度に効いている** (2b で言う「SN で遅くなる」の実体)。
  dt の戻りが 1 step 1.1 倍なのはコード側の制限ではなく (`newdt_fine`/`courant_fine` にそういう制限は無い)、
  熱いセルが膨張して音速が下がる物理的な時間。対策候補: SN の熱エネルギーを入れるセルの最小質量
  (軽いセルには隣接セルにも配る) か最高温度の上限。物理を変えるので本番前に要検討。

### 2g. 解法の自動選択 `cg_auto` (2026-11-11 追加)

grid 数だけでは CG/MG のどちらが速いか決まらない (CG の反復数は refine 領域のセル単位の差し渡しで決まり、
同期律速かどうかはランク数で決まる) ので、走りながら測って選ぶ。`patch/takashi/poisson_auto.f90`。

- `&POISSON_PARAMS cg_auto=.true.` (既定 .false.)、`cg_auto_every=20`、`cg_auto_margin=0.1`。
  `.false.` なら従来どおり `cg_levelmin` で固定 (旧 amr_step の分岐と bit 一致)。
- `.true.` のとき `cg_levelmin` は初期値。level ごとに「試行フラグ」`cg_trial_due` を持ち、ランの開始時と
  `cg_auto_every` coarse step ごとに全 level に立てる。フラグが立った level は次の解法を「今選んでいない方」で
  行い、同時に選んでいる方の計測値を捨てる。その次の解法 (普通は同じ coarse step の次の fine step) で
  選んでいる方を計り直し、両方揃った時点で MPI_ALLREDUCE(MIN) した値を比べ、`cg_auto_margin` 以上速ければ
  切り替える (ヒステリシス)。MIN を取るのは、最後に到着したランクの時間が待ちを含まない純粋な解法時間だから。
  決定は全ランクで同じ reduce 後の値で行うので分岐は一致する。
  rank 1 が ` cg_auto: level 13  MG 0.0123 s  CG 0.0080 s  -> CG` を出す。
- **まだ存在しない level** (宇宙論ランの初期など): `amr_step` は grid の無い level を呼ばないので、その level の
  試行フラグは立ったまま残り、level が初めて現れた最初の解法で試行・判定される。それまでは `cg_levelmin` の
  初期値 (既定 999 = MG) が使われるだけで害はない。現れたばかりで grid が少ない level の判定は意味が薄いが、
  次の試行 step で見直される。
- CG が `itermax=10000` で打ち切られたら (`phi_fine_cg.f90` に `cg_failed(ilevel)=.true.` を追加) その level は
  以後 MG 固定 (ログに 1 行)。
- 切り替えると解が丸め誤差レベルで変わる (どちらも epsilon まで収束した解だが残差のノルムが違う)。
  bit 再現性が要るランでは `.false.`。
- 状態 (`use_cg`, `t_solve_mg/cg`, `cg_failed` など) は `poisson_commons` に置き `init_poisson` で確保。restart では
  初期化し直されるので最初の coarse step で再判定する。モジュール変更なので `bin/*.o *.mod` を消してビルド。
- ファイル: 新規 `poisson_auto.f90`、patch に新規コピー `poisson_parameters.f90`, `phi_fine_cg.f90`、変更
  `poisson_commons.f90`, `init_poisson.f90`, `amr_step.f90`, `read_params.f90` (to_dice も), `bin/Makefile*` (POISSONOBJ)。
- **本番テスト flat-339614** (896 ランク、`cg_auto=.true.`, `cg_auto_every=1`, cg_levelmin 既定、2026-11-11):
  最初の coarse step で level 11–15 が CG、8–10 が MG に決まり、3 step とも同じ判定 (解法 1 回の時間は
  L15 MG 0.080–0.087 s / CG 0.039–0.041 s、L12 0.061–0.063 / 0.040–0.046、L10 0.057–0.059 / 0.058–0.064 で拮抗、
  L8 0.107 / 0.138)。poisson は 1 step あたり 34.2 / 30.5 / 29.8 s で cg_levelmin=13/11 の手動設定 (338146/338148) と同じ。
  判定の行は coarse step の先頭ではなく各 level の 2 回目の解法の直後 (level l なら 2^(15−l)+1 番目の fine step) に出る。
- **バグ修正 (同日)**: `nsubcycle(levelmin)=1` だと level levelmin+1 は 1 coarse step に 1 回しか解かれず、
  `cg_auto_every=1` では試行 (CG) → 次の step でまた試行 → … と選んでいる方 (MG) を一度も計らず、level 7 が
  毎 step CG (131 反復) のままになっていた (339614 で level 7 の判定行が無い)。判定待ち (`cg_pending`) の level には
  試行フラグを立て直さないようにした: 試行の次の解法で必ず選んでいる方を計って判定する。
  `cg_auto_every>=2` では起きない。

## 3. 通信の隠蔽 (Phase 2) のやり残し

- 実装済み: 全変数を 1 メッセージに集約 (`make_virtual_fine_dp_nv`)、`unew`/`rtunew`
  の reverse 通信を godunov と重ねる非同期化 (`make_virtual_reverse_dp_start/end`)。
  forward 版の非同期化は `make_virtual_fine_dp_start/end(xx,nvtot,iv1,iv2,ilevel,ih)` (スロット `async_fwd`,
  tag 210+ih)。start 時点の値を送信バッファに詰めるので start 後に xx を書き換えてもよいが、
  end までは仮想セルが古い値のまま。
- 実装済み (2026-10-30, bit 一致): `star_formation` 冒頭の `uold(1:4)` のブロッキング交換も同じ仕組みに。
  `rt_step` の最後の subcycle で cooling part 1 の直後に slot 1 で start (RT を使わない経路では粒子移動の
  直前に start)、`star_formation` が end してから ghost セルを active セルと同じ式 (`u=uold(2)/d`) で
  速度に変換する (元は primitive に変換してから送っていた。同じ数に同じ演算なので bit 一致)。
  通信が cooling part 2 と `move_fine` に隠れる。手元では star formation 0.21 s のまま。
- 実装済み (2026-10-30, bit 一致): RT subcycle 内の `rtuold` の forward 通信を cooling と重ねる。
  `cooling_fine_part(ilevel,1)` が他ランクに送る grid (emission grid) だけ先に冷却し、
  `make_virtual_fine_dp_start(rtuold,...,2)` で送信を出してから `cooling_fine_part(ilevel,2)` で
  残りを冷却、`make_virtual_fine_dp_end` で受信を仮想セルに書く (`amr_step.f90` の `rt_step`)。
  `coolfine1` はセル単位で独立なので grid の順序を変えても出力は bit 一致 (手元 8 ランク、96 ファイル)。
  手元 (1 ノード) では pack/unpack の帯域が支配的で時間は変わらない。本番 (rt - ghostzones 148 s/3 step,
  min 65 s) で純粋な通信分がどれだけ隠れるかを次の timer 表で見る。スロット: `async_fwd(2)`。
- 未実装: step 末の `uold` の forward 通信 (ghost 更新) の隠蔽。直後の `flag_fine`
  (勾配判定) が ghost を必要とするため、cooling_fine を「送信対象 grid → それ以外」の
  順に分けて送信を先に出す、または flag_fine を interior/boundary に分ける必要がある。
  star_formation の順序も絡むので物理順序の確認が必要。
- 未集約の per-variable ループ: `adaptive_loop.f90` (初期化時), `init_refine.f90`,
  `poisson/multigrid_*` (MG 独自の通信), `rho_fine.f90`, `pstarold/rho_eq/p_eq` など
  単発の変数。頻度が低いので優先度低。
- `make_virtual_reverse_dp` の `switchlevel=3` 以下 (粗い level) の rendezvous 方式は
  非同期版では使っていない (常に IRECV/ISEND)。超多ランクで粗い level のメモリが
  問題になれば検討。

## 4. 運用上のメモ

- **並列度を変えた restart (2026-10-29 確認)**: RAMSES 自身は書いたときと同じ ncpu でしか restart できないが、
  Mera.jl (v1.10) の `repartition` が output を別の ncpu 用に書き直せる。手元の 8 ランク output_00002 を
  4 ランクに書き直し、`compareramses` で octree・hydro/grav/rt の値・粒子 (id, x, v) が完全一致、
  `mp0`/`msp0`/`tp`/`zp`/family/tag も全粒子で一致、4 ランクで restart して完走した。
  手順 (Julia 1.13 は `/opt/homebrew/bin/julia`、`/usr/local/bin/julia` は古い 1.5 のまま):
  ```
  julia --project=<env> -e 'using Pkg; Pkg.add("Mera")'
  julia --project=<env> -e 'using Mera; r=repartition(51, "/path/run", 512; outpath="/path/run_512"); println(r);
                            println(compareramses(51, "/path/run", "/path/run_512"))'
  mpirun -np 512 ramses3d run.nml   # nrestart=51, ngridmax/npartmax は r の ngridmax_needed/npart_max 以上に
  ```
  注意: `stars_*.out` (SN ログ) はランク別ファイルなのでコピーされない (警告が出るだけで実害なし)。
  乱数列はランク数で変わるので星形成・SN は統計的な同等まで。QUADHILBERT ビルドは非対応。
  `part_file_descriptor.txt` で `mp0` と `msp0` が両方 `initial_stellar_mass` だったのを
  `stellar_population_mass` に改名した (`output_part.f90`、バイナリは不変)。名前で読むツール向け。
- `patch/takashi/` には未追跡ファイルが残っている (`rt_cooling_module.f90.back`,
  `to_dice/init_flow_fine_corrected.f90`, `backup_*` ディレクトリ, `format_fortran_script*.sh`)。
  必要なものは `git add`、不要なものは削除して状態をはっきりさせたい。
- 2026-10-29: per-level 分割のために patch に置いていた `load_balance.f90, refine_utils.f90,
  sink_particle.f90, flag_formation_sites.f90, init_tracer.f90` は削除した (上流が使われる)。
