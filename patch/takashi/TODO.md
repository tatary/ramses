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
