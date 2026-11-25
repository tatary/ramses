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

### 1b. (修正済 2026-11-19) `to_dice/init_part.f90`: `zp` の確保条件
- `zp` (粒子金属量) は `if(star.or.sink)` の中で `if(metal)` のときだけ確保されていたが、
  `to_dice/init_flow_fine.f90` は `metal` だけ見て `zp` を読む。`eps_star=0` で `metal=.true.` だと
  未確保配列を読んで SIGSEGV (変更前バイナリで再現)。
- 修正は 3 か所 (どれも `to_dice/` だけ): (1) `init_part.f90` の `zp` の確保を `star.or.sink` の外へ、
  (2) `load_dice` の `zp(ipart)=zz(i)` を `if(star)` の外へ、(3) `particle_tree.f90` で
  `dice_init.and.metal.and..not.(star.or.sink)` のとき `zp` をランク間で送る (幅 +1; 星ありのときは従来どおり
  `tp` の次で送るので変わらない)。(3) が無いと落ちはしないが、ランクを移ったガス粒子の金属量が 0 になる。
- 確認 (DICE テスト 4 ランク 6 step): 星ありは変更前と bit 一致。`eps_star=0` は完走し、初期の金属量の
  総和・min・max が星ありと一致 (1 ランクとも丸めの範囲で一致)。`LIGHT_MPI_COMM=1` でも heavy と bit 一致。

### 1c. (修正済 2026-11-19) RT ビルドで `rt=.false.` かつ `neq_chem=.false.`
- `rt_init.f90` の `read_rt_groups` で ion の添字 (`ixHII`) が 0 のまま `group_csn(…,0)` / `group_cse(…,0)` に
  **書き込んで**いた (範囲外)。bounds check 付きでは落ちる。通常のビルド (GNU) では落ちずに完走し、
  平衡冷却 (`cooling_module`) で進む。
- 修正: 該当の 6 行に `if(ixHII .gt. 0)` を付けた。`rt=.false.`、`neq_chem=.false.` のランは変更前と bit 一致
  (範囲外の書き込みは結果に効いていなかった)、既定テストも bit 一致。
- 一度 `neq_chem` を `.true.` に強制する形で入れたが (aebc95fe)、RT ビルドで平衡冷却が選べなくなるので外した。

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

### 1g. (修正済 2026-11-12) `rt_cooling_module.f90` の self-shielding が 3 < z < 10 で NaN、cloudy テーブルが z > 10 で範囲外

- `ss_factor_rahmati` (cloudy_metal_cooling かつ self_shielding のときの UVB 遮蔽) が α1, α2, β (負の値) を
  `linear_interpolation_in_log` で補間していて `log10(負)` = NaN。Chardin の格子 (z=3〜10) の間の z で必ず起き、
  UVrates 経由で化学・冷却全体に広がる。z ≤ 3 (DICE の z=0 run) と z ≥ 10 は格子端の値をそのまま使うので無事。
  テーブルを作った `ramses_tables/cloudy_cooling/maincl.cpp` と同じく α, β は線形補間、n0 と f だけ log 補間に
  変えた (`linear_interpolation` を追加)。z ≤ 3 は bit 不変。
- `rt_cmp_cloudy_metals` は `z ≥ z_reion` で CIE テーブルに切り替えるが、UVB 付きテーブルは z ≤ 10 (`zred_uv_max`)
  までしかない。`z_reion > 10` を指定すると 10 < z < z_reion で `izred ≥ 21` となり配列外を読む。判定を
  `z ≥ min(z_reion, zred_uv_max)` にした。既定の `z_reion=8.5` では不変。
- テーブル側の確認: `maincl.cpp` は Fortran と同じ Chardin et al. 2018 の z=3〜10 の値で、z < 3 は z=3 に固定、
  z=0〜10 を log(1+z) 等間隔 20 bin + CIE 1 bin (`zred=-1`)。`dlogz` の計算も `generate_fortran_tables.py` と一致。
  同じディレクトリの `cloudy_cooling.cpp` は別の格子・別のパラメータ (z<3 に Rahmati 2013 元論文の値) を使う古い
  生成コードで、今のテーブルには使われていない。
- UVrates の z_reion 付近の減衰 (`inp_UV_rates_table(..., z_damp=.true.)`、`10^(-20 (z/z_reion)^6)` 倍) は
  cloudy テーブル (z < z_reion は減衰なしの HM12) と合わないので `.false.` にした。`inp_UV_groups_table` 側は
  以前から `.false.`。これで H/He も金属も z_reion で一斉に UVB が入る。z=0 では減衰係数が 1 なので DICE run は不変。

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

- **2 回目 (flat-339648、`--hint=nomultithread`、1 ランク物理コア 4 個、`OMP_NUM_THREADS=4`)**: Intel MPI の pin 表は
  rank 0 が {0-3,112-115}、rank 1 が {56-59,168-171} で意図どおり 4 コア。しかし結果は 339635 (2 コア + HT) と
  ほぼ同じ (465 s 対 473 s、hydro-godunov 23.4 対 25.3 s/step、rt-godunov と cooling は同一)。
  **原因は binary に OpenMP が入っていなかったこと**: `Makefile.xd_intel` は `COMPILER = XD-INTEL` なのに、
  `OPENMP=1` の `-qopenmp` を `GNU`/`INTEL` の分岐にしか足していなかった (2026-11-11 に `XD-INTEL`/`XD-GNU`/`XD-CRAY`
  の分岐にも追加、`Makefile.xd`/`Makefile.xd_dice_intel` にも `OPENMP`/`BITREPRO`/`COOLRATES_NOCACHE` を追加)。
  つまり 339635 と 339648 はどちらも「束縛を明示した flat」で、339614 に対する 1.4 倍は **束縛だけの効果**。
  OpenMP の正味の効果は、`-qopenmp` 付きで作り直した binary (`ldd ramses3d | grep -i iomp` で libiomp5 が
  出ることを確認) で `NumThreadsPerProc=4` と `=1` を流し直して測る。
- **flat 運用への示唆**: 束縛を明示する (`--hint=nomultithread`、`--cpu-bind=cores`) だけで 3 step が 674 → 465 s
  (1.45 倍)。`run_partial_intel.rb`/`run.rb` にも `--hint=nomultithread` を入れた (b83418e9, fd825fa3)。
  339614 で何が起きていたか (HT の兄弟に 2 ランク?) は `--cpu-bind=verbose` を付けて流せば分かる。

- **3 回目 (flat-339656、`-qopenmp` 入り binary、`OMP_NUM_THREADS=4`、4 コア/ランク、2026-11-11)**: 今度は効いた。
  3 step 合計 465 → 296 s (束縛済み flat 339648 比 1.57 倍、元の 339614 比 2.28 倍、10/29 の 1258 s 比 4.25 倍)。
  coarse step は 115 / 90 / 88 s。1 step あたり (339648 → 339656):

  | 部分 | 339648 (1 thr) | 339656 (4 thr) | 比 |
  |---|---|---|---|
  | hydro - godunov | 23.4 s | 9.3 s | 2.52 (3 step とも) |
  | rt - godunov (step 3, subcycle 256 同士) | 8.0 | 2.8 | 2.8 |
  | cooling (step 2) | 16.1 | 5.8 | 2.8 |
  | hydro - rev ghostzones | 16.3 | 7.1 | 2.3 (待ちが減った) |
  | poisson | 17.6 | 23.6 | **0.75 (遅くなった)** |
  | rho / refine | 2.2 / 4.5 | 2.8 / 5.0 | 0.8 / 0.9 |
  | star formation | 17.6 | 17.5 | 1.0 |

  RT subcycle と dt 崩壊の履歴は 2 ランで違う (dt<1e-5 の fine step: 339648 は 0/2/246、339656 は 157/215/0) ので
  rt/cooling は履歴の揃った step で比べた。hydro-godunov は履歴に依らず 2.52 倍 (手元 Mac の 3.2 倍より低い:
  28 ランク × 4 スレッドでメモリ帯域を共有)。
- **スレッド化していない部分が 10–35% 遅くなった** (poisson 17.6 → 23.6 s、rho、refine)。解法 1 回の時間でも
  MG が L15 0.053 → 0.084 s (CG は 0.023 で同じ)。原因の候補は Intel OpenMP の既定 `KMP_BLOCKTIME=200` (ms):
  parallel region が終わった後もワーカーが 200 ms スピンするので、master が逐次部分を走る間も 112 コア全部が
  busy のままになり、all-core turbo で周波数が下がる (1 コアだけ busy なら高い)。試す: `KMP_BLOCKTIME=0`
  (`run_hybrid_intel.rb` に `KmpBlockTime` を追加、既定 "0"; "200" で元に戻る)。parallel region ごとの起床コストは
  数十 µs × 1 fine step あたり数回なので無視できる。
- **次の律速** (step 3、88 s): poisson 23.6、star formation 17.5 (min 11.8 / max 22.2: `uold` 交換の待ちを含む)、
  hydro-godunov 9.3、hydro-rev ghostzones 7.1、cooling 5.8、refine 5.0。poisson の multigrid smoothing / CG の
  matvec と、star formation の中身のスレッド化が次の候補。

- **OpenMP binary の 1 スレッド (flat-339657、339648 と同じ履歴)**: hydro-godunov 23.9 s (339648 の 23.4 と同じ)、
  rt-godunov、poisson も同じ → `-qopenmp` 自体の損は無い。**ただし cooling だけ 2 倍遅い**
  (20.7/16.1/8.1 → 42.6/33.1/16.7 s)。その待ちで rt-ghostzones と star formation も増えている。
  4 スレッドの cooling 5.8 s は「1 スレッド OpenMP 版 16.7 s の 2.9 倍」だが「非 OpenMP 版 8.1 s の 1.4 倍」でしかない。
  cooling の経路にだけある OpenMP 構文は `coolrates_module` の `inp_coolrates_table` の `threadprivate` キャッシュ
  (セルごと・反復ごとに何十回も呼ぶ最内側)、統計カウンタの `atomic`、初回の `critical`。量的に合うのは
  threadprivate (Mac の emulated TLS で見たのと同じ症状)。ifx の既定は `-qopenmp-threadprivate=compat` (native TLS)
  のはずだが、実際に遅い。**切り分け**: `Makefile.xd_intel` で `COOLRATES_NOCACHE = 1` にして
  `bin/coolrates_module.o` だけ消して再ビルド (他は変わらない)、1 スレッドで流して cooling が 339648 の 1.2 倍程度に
  戻れば threadprivate が原因 → 本番は NOCACHE=1 で行く (Mac では +17% の計算、2 倍よりずっと安い)。
  戻らなければ `atomic` か `-qopenmp` が暗黙に付ける `-auto` (全ローカルをスタックへ) を疑う。

- **`KMP_BLOCKTIME=0` (flat-339670、4 スレッド、NOCACHE=0)**: 狙いどおり逐次部分が戻った。poisson 23.6 → 18.2 s/step
  (339648 の 17.6 とほぼ同じ)、rho 2.8 → 2.3、refine 5.0 → 4.4、flag 2.3 → 1.7。代償はスレッドの起床:
  hydro-godunov 9.3 → 9.7 (+4%)、rt-godunov は subcycle 数の揃う step 1 (656) で 339648 比 2.66 倍 (BT=200 の
  step 3 では 2.9 倍)。差し引き 1 step あたり約 4 s の得なので **既定は 0 のまま**。中間の `KMP_BLOCKTIME=1` (ms;
  RT subcycle 内の連続する parallel region はカバーし、poisson の数十 ms の逐次部分では眠る) を試す価値はある。
  3 ラン (339656/339670) は dt 崩壊の履歴が違う (dt<1e-5: 157/215/0 対 0/0/37) ので合計時間 (296 対 310 s) では
  比べられない。履歴の揃う step 1 同士 (subcycle 656、崩壊なし) で: 339614 283.9 s → 339648 198.7 s (束縛) →
  339670 115.5 s (4 スレッド、BT=0): **元の flat 比 2.46 倍**。
  cooling は step 1/2 で 339648 比 1.35 倍しか伸びていない (threadprivate の件、NOCACHE=1 の切り分け待ち)。

- **`COOLRATES_NOCACHE=1` の 1 スレッド (flat-339672、339648 と同じ履歴)**: cooling 25.0 / 19.5 / 9.8 s で、
  非 OpenMP 版 (20.7 / 16.1 / 8.1) の 1.21 倍 (Mac で見た +17% と同じ)、threadprivate 版 (42.6 / 33.1 / 16.7) の半分。
  **cooling の 2 倍は ifx の `threadprivate` が原因と確定**。他の行は 339648 と同じ。
  **修正 (同日)**: キャッシュを `threadprivate` ではなく呼び出し側が持つ `type(coolrates_lookup)` にして
  `inp_coolrates_table(..., lookup=lk)` の optional 引数で渡す形にした。キャッシュが効くのは同じ T での連続した
  表引き (1 回の `compCoolrate` の中の 14 回) なので、`compCoolrate` のローカル `lk` で十分。`rt_cooling_module` の
  alpha の表引き (2〜6 回) は引数無し = 再計算のまま。TLS も save も無くなったので `COOLRATES_NOCACHE` は不要になり
  Makefile から削除。結果は旧キャッシュ版・NOCACHE 版と算術が同じなので bit 一致する (Mac、BITREPRO=1 の直列
  ビルドで参照と 40 ファイル bit 一致を確認)。Mac (gfortran、1 ランク、P コア) の cooling: 1 スレッド 25.4 s
  (旧 threadprivate 版 25.3、NOCACHE 版 29.6)、4 スレッド 7.34 s (3.46 倍; NOCACHE 版 7.98、threadprivate 版 12.85)。
  TOTAL は 1 スレッド 61.9 s → 4 スレッド 27.8 s。クラスタでは `coolrates_module` のモジュール I/F が変わるので
  `bin/*.o *.mod` を消してビルドし、1 スレッドで cooling が 339648 (20.7/16.1/8.1 s) に戻ることを確認する。

- **新キャッシュ版の本番テスト (2026-11-12)**: 1 スレッド flat-339730 (339648 と同じ履歴) は cooling 22.8 / 17.8 / 9.0 s で
  非 OpenMP 版の 1.10 倍 (threadprivate 版は 2.06 倍、NOCACHE 版は 1.21 倍)。残り 10% は `atomic` のカウンタと
  `-qopenmp` の `-auto` と思われるが追わない。他の行は 339648 と同じ。
  4 スレッド flat-339731 (`KMP_BLOCKTIME=0`、履歴 656/512/512、dt 崩壊なし): coarse step 104.0 / 90.9 / 88.8 s。
  履歴の揃う step 1 で 339648 比: hydro-godunov 23.4 → 9.6 s (2.43 倍)、rt-godunov 19.7 → 7.6 (2.58)、
  cooling 20.7 → 8.9 (2.32)、hydro-rev ghostzones 19.2 → 7.1、star formation 43.4 → 17.6、poisson 20.9 → 21.2 (同じ)。
  **step 1 の壁時計: 元の flat 339614 283.9 s → 束縛 339648 198.7 s → 4 スレッド 339731 104.0 s (2.73 倍)**。
- **次の課題: ランク間のばらつきが広がった**。step 1 の hydro-godunov が min 7.9 / max 14.3 s (1.8 倍; flat では
  20.8 / 26.4 = 1.27 倍)、cooling 7.7 / 13.8、rt-godunov 6.3 / 12.3。スレッドの効きがランクごとに違う
  (grid の多いランクほど parallel region が長く効率が良い) ので、leaf セル数 × 粒子数のコストモデルがずれる。
  待ちは rev ghostzones / rt-ghostzones / star formation に出る (star formation は min 11.6 / max 24.6)。
  **ばらつきの正体 (339731 step 2 の解析)**: hydro-godunov のランク間の差は grid 数と相関しない (corr 0.18)
  が、同じランクが cooling (0.83)、rt-godunov (0.91) でも遅い。ノード間の差ではなくノード内 (分散の 88%) で、
  **ノード内の task 位置 6, 7, 20, 21 が一貫して 10–17% 遅い** (位置別平均 9.7 s に対し 10.8 / 11.4 / 10.8 / 10.7;
  1 スレッドの 339648 では差なし)。この 4 つは物理コア 12–15、40–43 (socket 0) と 68–71、96–99 (socket 1) =
  **56 コアを 14 コアずつ 4 分割した境界をまたぐ 4 コアブロック** (Sapphire Rapids の tile / SNC4 の NUMA 境界)。
  28 task × 4 コアでは 14 が 4 で割れないので必ず 4 task がまたぐ。
  **`numactl -H` で確認済み (2026-11-12)**: NUMA node 8 つ (SNC4)、node0 = コア 0–13 (+HT 112–125)、node1 = 14–27、
  … node7 = 98–111。距離は同一 node 10、socket 内 14、socket 間 23。位置 6/20/7/21 の 4 コアブロック
  (12–15, 40–43, 68–71, 96–99) がちょうど node 境界をまたぐ。
  対策: (a) `cost_time_weighting=.true.` — 遅さがランクに固定された性質なので実測時間の補正で吸収できるはず
  (2c の仕組み; namelist だけで試せる)。(b) NUMA に揃えたランク配置: 14 コア = 3 task × 4 コア + 2 コア空き
  (24 task/ノード = 768 ランク) か 2 task × 7 スレッド (16 task/ノード = 512 ランク)。ランク数が変わるので
  Mera.jl の repartition が要る。(c) またぐ 4 task だけスレッド数を減らす wrapper (SLURM_LOCALID で分岐) は hack。
  step 1 の内訳 (104 s): poisson 21、hydro-godunov 9.6、cooling 8.9、rt-ghostzones 9.7、rt-godunov 7.6、
  star formation 17.6、hydro-rev gz 7.1、refine 5.1、rt-rev gz 5.0。
- **`cost_time_weighting=.true.` の 4 スレッド (flat-339758、2026-11-12)**: 仕組みは設計どおり動いた
  (time cost factor min 0.87 / max 1.36、前 step の実測効率と corr 1.00、またぐ 4 位置の grids-work が 0.90–0.93
  に減り、その hydro-godunov は 10.8/11.4/10.8/10.7 → 10.2/10.5/10.0/10.0 s) が、**全体は良くならない**:
  履歴の揃う step 3 (subcycle 512、崩壊なし) で 88.8 → 93.0 s (+4.7%)。3 カーネル合計のランク最大/平均は
  1.34 → 1.37 で変わらず (factor>1.15 の 56 ランクは仕事 0.88 倍でも効率 1.19 倍なので合計はまだ平均の 1.05 倍)、
  一方 grid 数の分割が 1.10 → 1.17 と不均一になったので grid 数に比例する逐次部分 (poisson +1.3 s、
  star formation +1.5 s、これらは同期を含む) が遅くなった。**`.false.` のまま** (2c の結論を維持)。
  改良するなら factor をスレッド化部分のコスト (全体の 1/4 程度) にだけ掛ける形だが、得られても数 s/step。
  ランクごとの効率差 (hydro 単独で 0.83–1.28、step 間で corr 0.82 と持続的) のうち NUMA またぎで説明できるのは
  corr 0.63 分で、残りは grid 数・level 数・ghost 比・星数のどれとも相関しない (1 スレッドでも 0.87–1.13 の持続的な
  差があり、4 スレッドとの corr は 0.49)。原因不明のまま保留。
- **`restart_remap=.true.` は step 1 を直さない (2026-11-12 に確認)**: 339614 以降は起動直後に "Load balancing"
  が出ている (restart_remap が効いている) のに、step 1 の分割は snapshot のものと grid 数で 0.2% しか違わず、
  grids-work max/avg は 1.69 のまま。step 2 の先頭の balance で 1.10 になる。起動時の balance は `init_tree` が
  粒子を levelmin までしか下ろしていない状態で走る (load_balance 自体は粒子を trunk に上げてから map を作るので
  粒子は直接は効かないはず) など、1 step 目の状態が in-run と違うことが原因と思われるが未解明。
  本番 (24 h = 数百 step) では 1 step 分なので優先度は低い。テストでは step 1 を比較に使わない。
  2026-11-14 の確認 (MUSIC zoom、16 ランク、単一分割の output から `two_tier_balance=.true., restart_remap=.true.` で再開):
  起動直後の balance で 2 層 (`TWO_TIER STATUS step 6 load balance ok`) になり、two-tier は最初の step から効く。
  起動時の分割と次の通常 balance の分割は粗い層でかなり違う (level 5 の active_max 935 → 2048) が、通常 balance 同士
  (step 7 → 8: 2048 → 1634、octs/cpu max 47752 → 34258) も同程度に違うので、restart_remap 特有の問題ではなく
  粒子支配の粗い層のコスト (1 oct + numbp×ngridmax/npartmax) に対する切れ目の敏感さ。起動時の粒子は in-run と同じく
  levelmin にある (init_tree は levelmin まで、in-run は merge_tree_fine で各 coarse step の頭に levelmin へ戻る)。

### 2d-4. poisson の OpenMP 化 (2026-11-15、`&RUN_PARAMS omp_poisson`、既定 .true.)

本番 (339731、4 スレッド、3 step 287 s) の poisson 59 s の内訳は log から CG (level 11-15、1 解法 50-80 反復、計 1500 回) が約 41 s、
multigrid (level 7-10) が約 3.4 s。level 15 は 1 ランク平均 2400 grid (1.9 万セル) で 1 反復 0.44 ms なので計算が主。

- `phi_fine_cg.f90`: 1 反復を 3 つの parallel region に (p の更新 / `cmp_Ap_cg` (orphaned `omp do`) と p·Ap / x と r の更新 + 次の r·r)。
  r·r は r の更新と同じループで求める (1 pass 減る; 値は同じ)。`cmp_Ap_cg` の save 作業配列を除去。
- `multigrid_fine_fine.f90`: 残差、ノルム、red-black Gauss-Seidel (同じ色のセルは独立)、restriction (1 つの fine セル番号では
  fine grid ごとに別の粗セル、8 個の加算順は逐次と同じ)、補間 (save 配列を除去)。**MG の核はスレッド数によらず同じ値**。
- `force_fine.f90`: `gradient_phi` のループと epot/rho_max のループ。`gradient_phi` と `interpol_phi` の save 配列を除去
  (`poisson/interpol_phi.f90` を patch にコピー; **未追跡なので commit 時に git add**)。
- 小さい level は逐次のまま (`omp_mg_mingrid=256` grid 未満、`amr_parameters`)。
- **和は `reduction` 句を使わない**: 部分和を合算する順序が実行ごとに変わり、CG では内積が次の反復に入るので同じスレッド数でも
  出力が再現しなかった (3 スレッド 2 回で 16/40 ファイル相違、丸め誤差の範囲だが回帰比較ができない)。スレッドごとの部分和を配列に入れて
  スレッド番号順に足す形 (`part(0:255)`) にした。CG の内積、MG のノルム、`force_fine` の epot (amr ファイルのヘッダに入る) が該当。
  hydro/RT の atomic (粗セル補正) は従来どおりスレッド数固定でも再現しない。
- 確認 (DICE levelmax 8、6 step、BITREPRO=1; CG は `cg_levelmin=7`): OPENMP=0 ビルドと OPENMP=1 の 1 スレッドは MG/CG とも変更前と 40 ファイル
  bit 一致。3 スレッド (omp_hydro/omp_rt=.false.): MG は omp_poisson の on/off で amr ファイル 4 個 (epot) 以外一致、CG は 2 回流して
  40 ファイル一致、1 スレッドとは総和・min/max の相対差 6.6e-10、反復回数は同じ。
- 手元 (1 ランク、コア固定、Xeon E5-2687W v4) の `poisson` timer: CG モード 15.6 → 6.5 s (4 スレッド、2.4 倍) → 5.5 s (7 スレッド、2.8 倍)、
  MG モード 12.5 → 6.7 s (1.9 倍) → 6.0 s (2.1 倍)。CG の解法だけなら 12.3 → 5.5 s (4 スレッド)。
- **本番で試すときは `KMP_BLOCKTIME=1`** (`run_hybrid_intel.rb` の `KmpBlockTime`): CG の region は 1 つ 0.1 ms 程度で、間に通信が入る。
  0 だと region ごとにワーカーを起こすコストが勝つはず (未測定)。1 ms なら CG の中では回り続け、長い逐次部分では眠る。
- poisson の残り (逐次): MG の粗い level (`multigrid_fine_coarse.f90`、上流のまま)、MG の準備 (`build_parent_comms_mg`、`make_fine_mask`、
  `make_fine_bc_rhs`、`restrict_mask`、`set_scan_flag`; 1 解法 4 反復なので比重が大きく、MG モードの伸びが 2 倍で止まる理由)、
  `make_initial_phi`、`cmp_residual_cg`、`save_phi_old`、`synchro_hydro_fine` (手元で poisson の 3%)。
- **未スレッド化の一覧 (339731 の割合)**: star formation 16.9%、rt-ghostzones 8.2%、hydro-rev ghostzones 7.2%、refine 5.0%、
  rt-rev ghostzones 3.9%、rho 2.5%、load balance 2.4%、flag 1.8%、feedback 1.4%、courant 1.2%、その他 3.5%。
- **star formation は計算ではなく待ちの可能性が高い**: 手元 (dice_rt_amr、4 ランク) では計算 0.06 s に対して `uold` 交換の完了待ち
  0.04-0.25 s、ALLREDUCE 5 本 0.001-0.07 s でランク間のばらつきが大きい。本番でもランク間で 35-66 s とばらつき、他の部分が
  スレッドで速くなると一緒に縮んだ (339648 → 339731 の step 1 で 43.4 → 17.6 s)。確かめるため `star_formation.f90` の timer を
  `star formation` (計算) / `sf - ghostzones` (交換の完了 + ghost の変換) / `sf - reduce` (集団通信) に分けた (出力は bit 一致)。
  次の本番ランの timer 表でこの 3 行を見る。待ちが大半なら star formation のスレッド化は不要で、負荷の不均衡の問題。

- **2026-11-17: クラスターで CG が遅くなっていた (896x4)**。`cg_auto` の 1 解法の時間 (level 15 / 14 / 13) は旧バイナリ f5ba86de
  0.0237 / 0.0296 / 0.0372 s、e14a460a の `omp_poisson=.false.` (340573) 0.0291 / 0.0358 / 0.0455 s (+22%)、スレッド化 + KBT 0
  (340628) 0.0419 / 0.0498 / 0.0541 s (+77%)。MG は 3 本とも同じ。タイマーの poisson は 59.2 / 68.3 / 76.4 s (3 step)。
  512x7 ではスレッド化で 86.4 → 65.4 s と効く (340461 → 340611)。
  スレッド化を切っても遅いので、疑いは反復ループの `collapse(2)` (flat index から (ind,i) を戻す形になる)。外して
  「各スレッドが ind ループを回し内側の i ループだけ分け合う (nowait)」形にした (`phi_fine_cg.f90`)。手元 (gfortran、
  dice_rt_amr + cg_levelmin=7、4 ランク): 1 スレッドは変更前と bit 一致で poisson 4.16 → 3.68 s、4 スレッドは 2.09 → 1.98 s
  (分け方が変わるので 4 スレッドの結果は丸めの範囲 (相対 1e-13) で変わる)。**ifx で戻るかは未確認**: クラスターで
  `cg_auto` の CG 時間が旧バイナリの値に戻るかを見る。896x4 の level 15 は 1 ランク 2400 grid で 1 反復 0.5 ms 弱なので、
  戻ってもスレッド化の得は小さい見込み (KBT 1 の結果待ち)。
- **KBT 1 の結果 (340569 `hu_lv15_omp4`、896x4、e14a460a、スレッド化あり)**: step 75.3 / 74.7 s、poisson 58.8 s、TOTAL 246.1 s
  (旧バイナリ 340011: 74.3 / 73.2 s、59.2 s、244.2 s)。CG 時間は level 15 / 14 / 13 / 12 = 0.0278 / 0.0338 / 0.0368 / 0.0265 s で、
  スレッド化が 1 スレッドの遅さ (+22%) をほぼ打ち消して旧バイナリ並み。KBT 0 (340628) の 76.4 s は KBT が原因と確定
  (並列領域が短いのでスレッドを寝かせてはいけない)。新バイナリの 896x4 の順位: KBT 1 + スレッド化 (75 s) < `omp_poisson=.false.`
  (78–79 s) < KBT 0 (81–82 s)。512x7 (340611) は 92.3 / 77.8 s、TOTAL 275.7 s。**本番は 896x4、KBT 1、`omp_poisson=.true.`**。
  collapse を外した版 (48206e6a) で 1 スレッドの遅さが消えれば、スレッド化の分だけ旧バイナリより速くなるはず (未確認)。
- **collapse を外した版 (48206e6a) の結果 (2026-11-18、896x4、KBT 1)**: 予想と逆で、**1 スレッドは直ったがスレッド化した CG が遅くなった**。

  | | step 2 / 3 | poisson | TOTAL | CG level 15 / 14 / 13 / 12 (s) |
  |---|---|---|---|---|
  | 旧バイナリ (340011) | 74.3 / 73.2 | 59.2 | 244.2 | 0.0237 (level 15) |
  | e14a460a `omp_poisson=.false.` (340573) | 79.0 / 77.9 | 68.3 | 258 | 0.0303 / 0.0391 / 0.0437 / 0.0325 |
  | e14a460a スレッド化 (340569) | 75.3 / 74.7 | 58.8 | 246.1 | 0.0281 / 0.0357 / 0.0377 / 0.0280 |
  | 48206e6a `omp_poisson=.false.` (341477) | **73.4 / 72.5** | **57.1** | **242.8** | 0.0237 / 0.0310 / 0.0343 / 0.0275 |
  | 48206e6a スレッド化 (341474) | 79.2 / 78.2 | 69.2 | 257.2 | 0.0411 / 0.0488 / 0.0515 / 0.0399 |

  1 スレッドの CG は旧バイナリと同じ 0.0237 s に戻った (遅さの原因は collapse(2) で確定)。スレッド化した CG は collapse 無しだと
  1.5 倍遅い (1 ランク 2400 grid を 4 スレッドで分けるループを 1 反復に 8×数回立てるので、ワークシェアの開始コストが勝つ)。
  どちらの形でも CG はスレッド化しない方が速い。一方 multigrid の粗い level (7–9) はスレッド化で 0.10 → 0.07 s と速く、
  340569 で CG が 0.028 s なのに poisson が 58.8 s だったことから CG 以外 (multigrid, force_fine) のスレッド化の得は数 s ある。
  そこで **`&RUN_PARAMS` に `omp_cg` (既定 .false.) を足した**: CG の反復は `omp_poisson .and. omp_cg` のときだけスレッド化し、
  `omp_poisson=.true.` は multigrid / force_fine だけに効く。**本番は 896x4、KBT 1、`omp_poisson=.true.`、`omp_cg` は既定のまま**
  (poisson 50 s 台前半の見込み、**クラスターで未確認**; 確認までは 48206e6a の `omp_poisson=.false.` が実測最速)。
  検証 (`tests/dice_rt_amr.nml` + `cg_levelmin=7`、4 ランク、BITREPRO、`omp_hydro=.false.`): `omp_cg=.true.` の 4 スレッドは 48206e6a の
  4 スレッドと bit 一致、既定 (`omp_cg=.false.`) は 4 スレッド = 1 スレッド = 48206e6a の 1 スレッドで bit 一致 (各 16 ファイル)。
  つまり `omp_cg` を切ると poisson の結果はスレッド数に依らない。
- **`omp_hydro=.true.` の 4 スレッドは実行ごとに結果が変わる (原因確定 2026-11-18、競合ではなく丸め)**: 上の検証中に、48206e6a 同士でも
  hydro/grav/part の 16 中 12 ファイルが不一致になった (`omp_hydro=.false.` か 1 スレッドなら一致)。原因は `godfine1` の粗セルへの
  flux 補正 (`unew/divu/enew(ind_buffer)` の `atomic` 6 箇所): 1 つの粗セルに最大 6 個の fine grid が足し、その順がスレッドの実行順で
  変わる (2d の最初に書いた既知の制限そのもの)。確認: (1) 2 回の 4 スレッドランの全セル統計 (`hydro_stats.py`) は相対 1e-13–1e-14 で一致、
  セル数も同じ。(2) 補正のブロック全体を `!$omp ordered` で囲んでループ順に固定した診断ビルド (リポジトリには入れていない) は
  4 スレッド 2 回 = 1 スレッド = 48206e6a の 1 スレッドで bit 一致 (16 ファイル) → 順序以外の原因 (競合) は無い。
  `ordered` は sweep を事実上直列にするので本番には使えない。再現性が要るなら補正を sweep ごとのバッファに貯めて並列ループの後に
  ループ順で足す形にする (RT の `rtunew` 補正 2 箇所も同じ)。**丸めの範囲なので受け入れる (2026-11-18 の判断)、実装しない**。4 スレッドの bit 比較は `omp_hydro=.false.` で行うこと
  (RT が光子を運ぶテストでは `omp_rt=.false.` も)。
- **flat 3584 (341473、`ngridtot=4e8`、`npartmax=400000`)**: メモリは足りて完走。`Used memory` 1.56 GB/ランク、step 98.5 / 96.1 s、
  TOTAL 350 s (poisson 77.1、rt - ghostzones 27.1、load balance 20.2 s)。896x4 (72.5 s) より 3 割遅く、flat に戻す理由は無い。

### 2d-4b. CG を速くする: warm start (`cg_warm_start`) と時間の内訳 (2026-11-18)

- 動機: クラスター (48206e6a、`omp_poisson=.false.`) の `cg_auto` の時間を反復数で割ると level 13–15 で 0.46 ms/反復とほぼ一定、
  1 ランクの grid が 1/4 の flat 3584 では 0.60 ms とむしろ遅い。計算量に比例していないので、スレッド化ではなく反復数と通信を減らす
  のが筋 (RAMSES-yOMP 論文 arXiv:2411.14631 も CG の得は PPCG による反復数 1/1.86 で、反復ごとの全体通信が最大のボトルネックと書く)。
- **`&POISSON_PARAMS` の `cg_warm_start` (既定 .false.、切れば変更前と bit 一致)**: `phi_fine_cg` の初期値を粗い level からの補間
  (`make_initial_phi`) ではなく前回の解にする (`make_warm_phi`)。
  - 前回の解は `phi_old` から取る。`phi` は解く直前に `rho_fine` が作業配列として上書きしている (最初 `phi` をそのまま使って
    反復が 3 倍に増えた)。`phi_old` は `amr_step` が `rho_fine` の直前に `save_phi_old` で写したもの。
  - 新しくできた grid には解が無い。`make_grid_fine` が 8 セルに親セルの値を入れているので、「8 セルが完全に同じ値」の grid だけ
    従来どおり補間値を使う (これで refinement 直後の反復が 38 → 27)。
  - 収束判定は warm のときだけ「残差 ≤ `epsilon` × 右辺ノルム」(従来は初期残差に対する比で、初期値が良いほど厳しくなる)。
  - 各 level のそのランで最初の solve は従来どおり (`cg_warm_ok`)。restart 直後も同じ。
- 結果 (`tests/dice_rt_amr.nml` + `cg_levelmin=7`、4 ランク、6 step): 反復数 level 7: 107 → 17、level 8: 101 → 25 (最初の 1 回を除く平均)、
  poisson 3.73 → 2.18 s。全セル統計は相対 1e-7–1e-6 の違い (ソルバーの許容誤差の範囲; スレッド化の丸めは 1e-13)。
- **未確認**: (1) 本番 (level 15、現状 52 反復) での減り方。(2) 許容誤差ぎりぎりの残差が同じ形で残り続けて長時間で系統誤差に
  ならないか (従来は毎回解き直すので誤差が入れ替わる)。ログの `econs` は冷却・フィードバック・RT 加熱を数えないので本番設定では
  判定に使えない; 冷却等を切った重力 + 流体のランか、従来の解との力の差で見る。(3) `cg_auto` との組み合わせ (CG が速くなると
  MG/CG の選択が変わる level が出るはず)。
- **CG の時間の内訳**: ランの最後 (タイマー表の直前、`clean_end` → `cg_timer_output`; `end.f90` を patch にコピー) に level ごとに
  solve 数、反復数/solve、set-up (初期値・右辺ノルム・最初の残差、solve あたり)、反復あたりの ghost 交換 / `MPI_ALLREDUCE` 2 回 /
  局所計算を、ランク平均 / 最大で出す。ALLREDUCE の時間には一番遅いランクを待つ分が入る。
  手元 (4 ランク) では局所計算が 0.7–1.1 ms/反復で支配的、**set-up が 20–35 ms/solve と 30 反復分**あり、warm start 後は solve の半分を
  占める (`interpol_phi` と `cmp_residual_cg`)。クラスターの内訳を見てから次の手 (set-up の軽量化、ALLREDUCE を 1 回にまとめる
  Chronopoulos–Gear 形、前処理) を決める。

- **set-up の軽量化 (2026-11-18)**: (1) `cmp_residual_cg` は束ごとに `interpol_phi` を 6 回無条件に呼んでいたが、補間値を使うのは
  その方向に隣の grid が無い grid だけなので、該当する grid だけ集めて呼ぶようにした (スレッド無関係、flat MPI にも効く)。
  (2) `cmp_residual_cg` / `make_initial_phi` / `make_warm_phi` の束のループを `omp_poisson` でスレッド化 (作業配列の `save` を外した;
  `make_initial_phi` は multigrid からも呼ばれる)。どちらも結果は bit 一致: 1 スレッドは変更前と (既定・warm とも 16 ファイル)、
  4 スレッドは 1 スレッドと (`omp_hydro=.false.`、既定・warm とも)。
  手元 (4 ランク、level 7 / 8、set-up ms/solve): 変更前 32.5 / 19.6 → (1) で 8.7 / 5.5 → 4 スレッドで 4.8 / 3.3。
  warm start では 28.2 / 17.1 → 4.7 / 2.8 → 2.7 / 1.9。CG の合計 (level 7+8) は 変更前の既定 2.23 s → warm + 軽量化 + 4 スレッド 0.57 s。

- **反復ループの形を solve ごとに分岐 (2026-11-18)**: スレッド化しないとき (既定の `omp_cg=.false.`、または grid が `omp_mg_mingrid` 未満) は
  OpenMP の指示を含まない素のループ、スレッド化するとき (`omp_poisson .and. omp_cg`) は collapse(2)。根拠はクラスター (ifx、896x4、
  level 15、1 ランク 2400 grid) の CG 時間: 1 スレッドは collapse なし 0.0237 s < あり 0.0303 s、4 スレッドは collapse あり 0.0281 s < なし 0.0411 s。
  検証 (4 ランク、16 ファイル): 既定の 1 / 4 スレッドは変更前と bit 一致、`omp_cg=.true.` の 1 スレッドは既定と bit 一致、4 スレッドは 2 回で一致
  (collapse なしの形とは和の分け方が変わるので丸めの範囲で不一致)。
  **手元 (gfortran、1 ランクの grid が多い) では 4 スレッドの collapse ありの方が 2 割遅い**: 反復の局所計算 level 7 / 8 = collapse なし 0.472 / 0.200 ms、
  あり 0.573 / 0.243 ms (1 スレッドは 1.02 / 0.63 ms)。クラスターと逆なので、コンパイラか grid 数のどちらに依るのかは未確認。
  `omp_cg=.true.` を本番で使う場面 (512x7 など grid が多い設定) が出たら、ifx でどちらの形が速いかを測り直すこと。
  なお手元では CG のスレッド化自体は 2.2–2.6 倍効いている (grid が多ければ `omp_cg` は得)。

- **warm start の精度 (2026-11-19)**: `tests/dice_rt_amr.nml` + `cg_levelmin=7`、4 ランク 1 スレッド、30 step、95ee9930。基準は `epsilon=1e-9` の
  従来の解き方。`grav_*` の力を同じセルどうしで比べ、差の rms を基準の力の rms で割った値 (level 7 / 8、5 step 目 = ランの軌道が分かれる前):

  | | 力の rms 誤差 L7 / L8 | phi の rms 誤差 L7 / L8 | 反復/solve L7 / L8 | CG 合計 (s) | poisson (s) |
  |---|---|---|---|---|---|
  | 従来 `epsilon=1e-4` | 2.6e-5 / 8.1e-6 | 2.4e-5 / 1.5e-5 | 107 / 102 | 9.8 | 15.9 |
  | warm `1e-4` | 3.6e-5 / 3.2e-5 | 2.5e-5 / 1.0e-5 | 22 / 30 | 2.6 | 8.5 |
  | 従来 `1e-5` | 3.3e-6 / 8.9e-7 | 3.4e-6 / 2.1e-6 | 144 / 136 | 12.5 | 19.0 |
  | warm `1e-5` | 4.4e-6 / 2.8e-6 | 3.2e-6 / 1.3e-6 | 50 / 70 | 5.6 | 12.1 |

  - 同じ `epsilon` だと warm は phi の誤差は同じか小さいが、**力の誤差は level 7 で 1.4 倍、level 8 で 3–4 倍大きい** (7 step 目も同じ比)。
    誤差が小さいスケールに残っている。それでも力の 3e-5 で、最大でも 2e-4。
  - **warm + `epsilon=1e-5` は従来の `1e-4` より全 level で正確 (力 6 倍 / 3 倍) で、CG の時間は 57%**。精度を落とさずに使うならこの組み合わせ。
  - 誤差の蓄積は見えない: 5 → 7 → 15 step で warm の rms 誤差は増えていない (L7 3.6e-5 → 3.6e-5 → 3.7e-5)。
  - 30 step 目はどのランも level 8 の力の差が 7–9e-5 に揃う (`epsilon=1e-5` の従来の解き方でも 7.7e-5)。ソルバーの誤差ではなく
    ランの軌道が分かれた分で、この長さより先は力の比較では判定できない。grid 構造は全出力で 5 本とも同じ。
  - run は `claude_scratch/run/a_{ref,cold4,warm4,cold5,warm5}`。

- **クラスター (95ee9930、896x4、KBT 1、3 step、`hu_lv15_omp4_nops/`、2026-11-19)**: A=342065 (`omp_poisson=.false.`)、B=342066 (A + `cg_warm_start`)、
  C=342067 (`omp_poisson=.true.`)、D=342219 (C + `omp_cg=.true.`)。比較は 341477 (48206e6a、A と同じ namelist)。

  | | 341477 | A | B | C | D |
  |---|---|---|---|---|---|
  | poisson (s) | 57.1 | 53.5 | 36.5 | 52.2 | 51.3 |
  | CG level 15: 反復/solve | – | 52.0 | 30.0 | 52.1 | 52.1 |
  | CG level 15: 合計 (s) | – | 17.5 | 10.6 | 18.3 | 17.9 |
  | CG level 7–15 合計 (s) | – | 40.6 | 23.8 | 41.6 | 41.2 |
  | level 15 の RT subcycle | 1011 | 1663 | 1540 | 1011 | 1011 |
  | TOTAL (s) | 242.8 | 280.6 | 253.4 | 239.9 | 238.4 |

  - **A と B は RT subcycle が 1.5–1.65 倍あり、TOTAL と step 時間は他と比べられない** (rt-godunov / cooling / rt-ghostzones / rt-rev / sf-reduce が
    その分増えている; 同じ namelist の 341477 と A で変わる理由は未調査)。この 5 項目を 341477 に揃えると A ≈ 242、B ≈ 223。
  - warm start: level 15 の反復 52 → 30 (level 10–14 は 115–68 → 34–42)、poisson −17 s (3 step)、全体の約 7%。手元 (107 → 22) ほどは減らない。
  - `omp_poisson=.true.` はもう遅くならない (set-up 2.64 → 1.95 ms/solve、poisson −1.3 s)。`omp_cg` の上乗せは −0.9 s で誤差の範囲
    (level 15 の局所計算 0.158 → 0.144 ms/反復、交換は 0.042 → 0.064 に増える)。
  - level 15 の 1 反復 0.39 ms の内訳 (A): ALLREDUCE 2 回 0.19、局所計算 0.16 (最大 0.22)、交換 0.04。set-up は 2.6 ms/solve で反復 7 回分。
    ALLREDUCE が最大なので、次の手は ALLREDUCE を 1 回にまとめる形 (Chronopoulos–Gear) か前処理で反復を減らすこと。
  - `epsilon` は MG の収束判定 (`multigrid_fine_commons.f90` の `err<epsilon`) にも使われる。本番 namelist は指定なし (既定 1e-4)。
  - 判定の分母: MG と従来の CG は補間の初期値の残差、warm は右辺ノルム。ログの `==> Level=` の 1 つ目 (残差 / 右辺ノルム) の最終値は
    level 15 で 従来 8.4e-5 / warm 9.4e-5、level 12–14 で 7.8–8.6e-5 / 9.6–9.7e-5、level 7–10 で 2.7–3.6e-5 / 6.7–9.2e-5。
    細かい level では初期残差 ≈ 右辺ノルムなので 2 つの基準はほぼ同じ、粗い level では warm が 3 倍ほど緩い。

- **判定の分母を初期残差にする案 (`cg_warm_ini_res`、2026-11-19 に試して外した; コードには無い)**: warm start の収束判定を、従来の CG・MG と同じ
  「補間の初期値の残差」に対する比にするオプション。solve ごとに補間 (`make_initial_phi`) と残差計算を 1 回余計に行って基準を求めてから、
  前回の解で始める。手元 (上の表と同じ 30 step のラン、`a_warmi4`、`epsilon=1e-4`): 力の rms 誤差 L7 / L8 = 2.7e-5 / 1.1e-5 で従来
  (2.6e-5 / 8.1e-6) に戻るが、反復 31 / 51 (warm 22 / 30)、set-up 12.9 / 8.0 ms/solve (warm 5.3 / 3.2; 全 grid の補間が戻る)、CG 合計 4.2 s (従来 9.8、warm 2.6)。
  - 基準を使い回す案 (初期残差 / 右辺ノルムの比を N solve ごとに更新) も検討した。比は連続する solve の間では 1% 未満しか変わらないが
    (クラスター level 15: 0.891–0.896)、比例する根拠はなく長時間ではずれる (手元 level 8: 30 step で 0.305 → 0.229)。
  - クラスターの反復は従来・warm とも残差 1 桁あたり約 13 回 (level 15)。warm を `epsilon` 3e-5 / 1e-5 にすると 30 → 約 37 / 44 回の見積もり (従来 52)。
  - **結論 (2026-11-19)**: 力の誤差 3e-5 は溜まらず (5 → 15 step で一定) 実害がないので、**本番は `cg_warm_start=.true.`、右辺ノルム基準、`epsilon` は既定 1e-4 のまま**。
    `cg_warm_ini_res` も `epsilon_cg` (CG だけ別の許容誤差) も入れない。`omp_poisson=.true.` と合わせた 1 本はまだ流していない。

- **zoom の宇宙論ランでは warm start は逆効果 (2026-11-19)**: `namelist/ks/l5_14_test.nml` (MUSIC zoom、levelmin 5 / levelmax 16、6 step、8 ランク) に
  `cg_levelmin=6` を足して従来と warm を比べた (`claude_scratch/run/zoom_zc`, `zoom_zw`; ビルドは `build_ks_new` = GNU、NGROUPS=5、NPSCAL=3、DICE 空)。
  反復/solve は level 6–8 で同じ (46–50)、level 10: 58 → 70、12: 86 → 130、13: 134 → 221、14: 176 → 375。CG 合計 9.2 → 17.2 s。発散や未収束は無い
  (最終 残差/右辺 は 7–10e-5)。
  - 前回の解の初期残差 / 右辺ノルムが level 7 で 3、level 10 で 63、level 14 で 1600 (従来の補間の初期値は全 level で 0.8–1.7)。
    level が 1 つ上がるごとに 2.5–3 倍になる。
  - このランは `nsubcycle` が level 5–13 で 1 (全 level が coarse step ごとに 1 solve)。subcycle のせいかを見るため level 10 以上を 2 にして
    2 step 流した (`zoom_zcs`, `zoom_zws`; level 14 は coarse step あたり 16 solve): それでも warm のほうが遅い。level 13: 134 → 179、14: 173 → 288 反復、
    CG 合計 49 → 68 s。solve ごとに見ると、最初の coarse step の中 (粗い level が動かない間) は前回の解の初期残差 / 右辺ノルム が 0.16–1.4 で
    反復は 172–257 (従来 173; 初期残差が補間より小さくても反復は減らない)、粗い level を解き直した直後は 600–5500 で約 400 反復。
    つまり solve の頻度 (subcycle) だけでは決まらない。
  - `cg_auto=.true.` (`cg_auto_every=1`、nsubcycle は元のまま 6 step; `zoom_zca`, `zoom_zwa`) だと level 6–12 は CG (MG 0.010–0.034 s に対し 0.005–0.028 s)、
    grid の多い level 13 / 14 は MG (0.086 / 0.55 s; CG は 0.117 / 1.11 s) が選ばれる。poisson は 9.3 s (`cg_levelmin=6` 固定の従来 10.8 s)。
    warm を足すと level 12 も MG になり、CG の試行が重くなって (level 14 で 1 回 2.95 s) poisson 15.0 s。8 ランクの手元の時間なので、選択はランク数で変わりうる。
  - 推定 (未検証): 高赤方偏移の zoom は密度がほぼ一様で右辺 (密度ゆらぎ × dx²) が小さい。level の端の境界値は粗い level の phi の時間補間で
    solve ごとに動き、その変化は phi の大きさで決まるので、右辺に対して大きな残差になる。これを消すのは level 全体にわたる滑らかな補正で、
    CG では反復がかかる (補間の初期値はその滑らかな部分を最初から持っている)。DICE は密度の濃淡が大きく右辺が大きいので、境界の変化が相対的に小さい
    (クラスター level 15 の初期残差 / 右辺ノルム は 0.015)。
  - 右辺ノルムが小さすぎて判定が成り立たない、という心配のほうは起きていない (ほぼ一様な level 6–9 でも初期残差 / 右辺ノルム ≈ 1)。
  - **当面: `cg_warm_start` は DICE の本番専用。宇宙論 / zoom では、少なくとも高赤方偏移では使わない** (低赤方偏移で銀河ができた後は未確認; 使う前に
    `==> Level=` の反復数を従来と比べること)。
  - 直すなら、初期値を「前回の解 + (今回の粗い level からの補間 − 前回の solve のときの補間)」にする。前回の補間値を覚える配列 (セルあたり実数 1 個)
    と、毎回の全 grid の補間 (set-up が従来並みに戻る) が要る。DICE でも反復がさらに減る可能性があるが未確認。

### 2d-4c. 次にクラスターで流すテスト (2026-11-19 時点; 実行済み、結果は 2d-4e)

バイナリは 95ee9930 のまま (作り直し不要)。場所は `synology1/dwarf/run/hu_lv15_heavy/hu_lv15_omp4_nops/`、`output_00051` から 3 step、
896 ランク × 4 スレッド、`KmpBlockTime "1"` (`run_hybrid_intel.rb` の最後の行で namelist を選ぶ)。namelist は `xd_rt_pos.nml` (`omp_poisson=.true.`、
`cg_auto=.true.`、`cg_auto_every=1`、`timer_ranks_every=1`) からの変更で書く。優先順:

| 優先 | 変更 | 目的 | 見るもの |
|---|---|---|---|
| 1 (必須) | `&POISSON_PARAMS` に `cg_warm_start=.true.` | 本番候補 (`omp_poisson=.true.` + warm)。この組み合わせは未実行 | poisson が 36 s 前後 (B=342066 は 36.5)、level 15 の反復が 30 前後、TOTAL |
| 2 | なし (C=342067 の再実行) | level 15 の RT subcycle 数が実行ごとに変わるか (A / B だけ 1663 / 1540、他は 1011) | `Performed level 15 RT-step` の subcycle の合計。1011 なら偶然、また増えるなら原因を調べる (RT と cooling が 1.6 倍になる) |
| 3 | 512 ランク × 7 スレッド (`hu_lv15_test_512/output_00051`、`CoresPerTask=7`、`NumThreadsPerProc=7`)、`cg_warm_start=.true.`、`omp_cg=.true.` | 4 スレッドと 7 スレッドの決着 (前回 340611 は 92.3 / 77.8 s で負け; grid が多いので `omp_cg` と warm が効くはず) | step 時間、poisson、CG の内訳。896x4 で本番を始めるなら省略可 |
| 4 | `cg_warm_start=.true.`、`cg_auto_every=20` | `cg_auto` の MG の試行 (3 step で 26 回) を減らす。得は 1 step 0.5 s 程度の見込み | poisson、`cg_auto:` の行数 |
| 5 | **新しいバイナリ** (`cg_precond` 入り; `.o` と `.mod` を消してビルド)、`cg_warm_start=.true.`、`cg_precond=1`、`cg_precond_omega=1.5` | 前処理 (2d-4d) の実測。他のテストとバイナリが変わるが、`cg_precond=0` なら 95ee9930 と同じ結果 | level 15 の反復 (30 → 14 前後の見込み)、CG の内訳の局所計算、poisson |

- どのログでも、TOTAL を比べる前に level 15 の RT subcycle の合計が揃っているかを確かめる。
- `epsilon` は既定 (1e-4) のまま。`omp_cg` は 3 以外では指定しない (既定 .false.)。
- **本番の namelist** は 1 の設定から `timer_ranks_every=1` を外す (既定 0)。`nstepmax` / `foutput` のテスト用の指定も戻す。
- リポジトリの `namelist/dice/xd_rt.nml` と `jobfiles/run_hybrid_intel.rb` は古い (`cg_auto` / `omp_poisson` なし、512x7、`KmpBlockTime "0"`)。本番設定が決まったら合わせる。

### 2d-4d. CG の前処理 `cg_precond` (2026-11-19、実験的、クラスター未測定)

- **`&POISSON_PARAMS` の `cg_precond=1` (既定 0、切れば変更前と bit 一致)**: 前処理付き CG。前処理は「自分のランクのセルだけで閉じた SSOR 1 往復」
  (`cg_precond_omega`、既定 1 = 対称 Gauss–Seidel)。他ランク・level の外の隣は 0 とみなすので通信は増えない (ランクごとのブロック Jacobi)。
  ALLREDUCE も 1 反復 2 回のまま (r·r と r·z を 1 回で送る)。収束判定は従来と同じ r のノルム。
  - 実装 (`phi_fine_cg.f90`): `cg_precond_table` が solve ごとに隣セルの表 (`cg_nb(6, セル)`、自分のセルでなければ添字 0 = 常に 0 の要素) を作り、
    `cg_precond_sgs` が前進 (`cg_yw`) → 後退 (`cg_zw`) の掃引をする。作業配列はセルあたり実数 2 個 + grid あたり整数 1 個 (最初に使うとき確保) と、
    その level の自分のセル × 整数 7 個の表。掃引は逐次で、`cg_precond>0` のとき反復ループはスレッド化しない形になる (`omp_cg` は無視)。
  - 結果は領域分割に依存する (ランク数や load balance で前処理が変わる; 許容誤差の範囲)。
- 手元 (`tests/dice_rt_amr.nml` + `cg_levelmin=7`、6 step、level 7 / 8、1 スレッド):

  | | 反復/solve | 局所計算 ms/反復 | CG 合計 (s) |
  |---|---|---|---|
  | 16 ランク 従来 | 107 / 101 | 0.265 / 0.164 | 0.73 |
  | 16 ランク 前処理 ω=1 | 53 / 49 | 0.464 / 0.286 | 0.63 |
  | 16 ランク 前処理 ω=1.5 | 52 / 42 | 0.463 / 0.286 | 0.58 |
  | 16 ランク warm | 30 / 31 | 0.267 / 0.165 | 0.24 |
  | 16 ランク warm + 前処理 ω=1.5 | 14 / 13 | 0.473 / 0.293 | 0.21 |
  | 4 ランク 従来 | 107 / 101 | 1.03 / 0.63 | 1.91 |
  | 4 ランク 前処理 ω=1.5 | 37 / 38 | 1.75 / 1.08 | 1.24 |

  16 ランクは level 8 で 1 ランク 1400 grid (本番 level 15 は 2300)。反復は 1/2.0–1/2.4 (4 ランクの大きいブロックでは 1/2.7–1/2.9)、
  局所計算は 1.74 倍。ω=1.7 は 1.5 より悪い (55 / 45)。表を作る前の版 (掃引のたびに `son(nbor)` を引く) は局所計算が 2.8 倍で損だった。
- **クラスターの見積もり (未測定)**: level 15 は 1 反復 0.39 ms のうち通信が 0.23 ms なので、局所計算が 1.74 倍でも 1 反復は 0.56 ms 程度、
  反復が 1/2.2 なら CG は 3 割強の短縮 (warm ありで 3 step 23.8 → 17 s 前後、全体の約 3%)。
- **力の精度 (2026-11-19)**: 2d-4b の精度検証と同じ 30 step のラン (4 ランク、`epsilon=1e-4`、基準 `a_ref` = `epsilon=1e-9`)。5 step 目の力の rms 誤差 L7 / L8:

  | | 力の rms 誤差 | 反復/solve | CG 合計 (s) |
  |---|---|---|---|
  | 従来 | 2.6e-5 / 8.1e-6 | 107 / 102 | 9.8 |
  | 前処理 ω=1 (`a_prec4w1`) | 2.6e-5 / 7.5e-6 | 49 / 47 | 7.8 |
  | 前処理 ω=1.5 (`a_prec4`) | 2.5e-5 / 5.5e-6 | 37 / 39 | 6.3 |
  | warm | 3.6e-5 / 3.2e-5 | 22 / 30 | 2.6 |
  | warm + 前処理 ω=1.5 (`a_wprec4`) | 2.9e-5 / 2.2e-5 | 8.5 / 12 | 2.0 |

  前処理は精度を落とさない (同じか少し良い; 細かいスケールの誤差を先に消すため)。7 / 15 step 目も同じ水準で、誤差は増えていない。
- **未確認**: (1) クラスターでの反復数と時間。(2) `cg_auto` との組み合わせ。(3) zoom (level 14 で 176 反復) での効き。(4) スレッド化 (スレッドごとのブロックにすれば掃引を並列にできる)。

### 2d-4e. 2d-4c のテスト結果 (2026-11-19 夜のジョブ、バイナリ 8511f0b9、2026-11-20 に整理)

ログは `synology1/dwarf/run/hu_lv15_heavy/hu_lv15_omp4_nops/` (896x4) と `hu_lv15_omp7/` (512x7)。level 15 の RT subcycle は 896x4 の 4 本とも
step 2 / 3 が 256 / 256 で揃っているので、step 2–3 の平均で比べる (`timer_ranks_*_n7576/8.dat` の差の半分、ランク平均)。

| テスト | ジョブ | 設定 | step 2 / 3 (s) | poisson (s/step) | level 15 反復/solve | CG 合計 (3 step, s) | RT subcycle (step 1/2/3) |
|---|---|---|---|---|---|---|---|
| 2 | 343084 | 変更なし (342067 の再実行) | 71.9 / 69.9 | 16.3 | 52.1 | 40.9 | 499 / 256 / 256 |
| 1 | 343082 | `cg_warm_start` | 65.2 / 64.8 | 10.3 | 26.9 | 23.8 | 499 / 256 / 256 |
| 4 | 343085 | warm + `cg_auto_every=20` | 63.8 / 61.5 | 7.9 | 19.5 (step 2 / 3 は 13.5 / 16.6) | 18.5 | 604 / 256 / 256 |
| 5 | 343086 | warm + `cg_precond=1`、ω=1.5 (`cg_auto_every=1`) | 63.6 / 62.2 | 8.6 | 11.1 | 18.2 | 495 / 256 / 256 |
| 3 | 343087 | 512x7、warm + `omp_cg` | 73.5 / 65.4 | 10.4 | 26.0 | 22.8 | 656 / 388 / 256 |

poisson 以外の行は 896x4 の 4 本で 0.1 s/step 以内で同じ (hydro - godunov 9.5–9.7、sf - reduce 8.3–8.4、hydro - rev ghostzones 6.5、refine 4.2、
star formation 3.6、cooling 3.5、rt - godunov 3.0、rt - ghostzones 3.1)。

- **warm start (1)**: 見込みどおり。level 15 の反復 52 → 27、poisson 16.3 → 10.3 s/step、step 時間 −8%。
- **RT subcycle のばらつき (2)**: 同じ設定の再実行は 342067 / 342219 と完全に同じ (499 / 256 / 256、反復数も同じ) なので実行ごとの乱れではない。
  設定を変えると (重力が許容誤差の範囲で変わると) step 1 の数が 495–656 と変わる。step 2 / 3 は 896x4 の今回の 4 本では揃ったが、
  342065 / 342066 は 495 / 512、416 / 504、512x7 は step 2 が 388。RT の dt を決めるセルの状態が丸め程度の違いに敏感ということで、
  3 step のテストでは TOTAL ではなく subcycle の揃った step 同士で比べるしかない。
- **512x7 (3)**: subcycle の揃った step 3 で 65.4 s、896x4 の warm (343082) は 64.8 s。poisson も 10.4 と 10.3 s/step で同じ。7 スレッドの得は無いので
  本番は 896x4 のまま。
- **`cg_auto_every=20` (4) は見込み (0.5 s/step) よりずっと効いた: poisson 10.3 → 7.9 s/step。** 理由は MG の試行の時間そのものではなく、
  **試行が warm start を乱すこと**。`cg_auto_every=1` (343082) では level 15 の反復 (64 solve ごとの平均) が coarse step のたびに
  45 → 26 → 24 → 16、40 → 35 → 17 → 12 と上がり直す。step 2 では level 9–15 のすべてが coarse step の前半ずっと 50–80 反復で、
  level 8 が 2 回目の solve (step の中間) をした直後に 10–25 に落ちる。試行の無い 343085 の step 2 / 3 でも step の最初の 1/8 (level 10 の 2 回目の
  solve まで) は 20–35 反復で、その後 6–20。
  - 推測 (未検証): 細かい level の境界値は粗い level の `phi_old` と `phi` の時間補間なので、粗い level の `phi − phi_old` に解法の違い
    (試行で MG ↔ CG を入れ替えた分の許容誤差程度の差) が入ると、その level が次に解かれるまで細かい level の境界が fine step ごとに動き続ける。
    1 step に 1–2 回しか解かない level 7 / 8 の試行が step の半分以上を乱す。試行の無い step の最初の 1/8 は、coarse step の頭で粗い level が
    解き直された分 (本物の時間変化) と思われる。
  - 343082 では step 2 の試行で level 8–11 が MG → CG に変わった (warm の CG は MG の解から 4 反復で収束するので試行だけ速く見える。
    切り替えた後は 50–70 反復)。343085 は level 7–10 が MG のまま。**warm と `cg_auto` の組み合わせでは試行の時間が実際の時間を表さない。**
  - 対策の候補: (a) 本番は `cg_auto_every` を既定の 20 以上にする (試行のある step だけ 343082 並み = +2.4 s)。(b) 試行の solve の後で選択中の解法で
    解き直して解を揃える (試行のたびに 1 solve 余計にかかるが、判定がその場で済み `cg_pending` も要らなくなる)。(c) CG の時間は試行の 1 回ではなく
    直近の数回の平均で比べる。手元 (`dice_rt_amr.nml` + `cg_levelmin=7`) で (b) を試す価値がある。
- **(b) を実装 (2026-11-20、`poisson_auto.f90`、クラスター未測定)**: 試行のある solve は「選択していない解法で 1 回 (時間を測るだけ)」→
  「選択中の解法でもう 1 回 (この解を残す)」と 2 回解き、その場で判定する。2 回目は同じ状態から始まる (初期値は `phi_old` か level-1 から作り、
  `rho` は solve で変わらない) ので、残る解は試行の無いランと同じ。`cg_pending` は削除 (`poisson_commons` が変わるので `.o` / `.mod` を消してビルド)。
  判定で解法が入れ替わった次の solve だけは従来どおり解が切り替わる。
  - 検証 (`run/a_warm4` の namelist を 10 step、4 ランク、level 7 / 8、`cg_auto_every=1`、80 ファイル): `cg_auto` + `cg_auto_margin=1` (切り替えなし) は
    `cg_auto=.false.` と bit 一致 (CG 選択 + MG 試行 = `rs_C` / `rs_D`、MG 選択 + CG 試行 = `rs_E` / `rs_F`)。既定の margin でも選択が変わらなければ
    一致 (`rs_B`)。変更前のバイナリは 48 ファイルが不一致で、CG の反復は 36.8 / 39.8 (変更後と試行なしは 25.5 / 29.4) と、乱れも手元で再現した。
  - 試行 1 回のコストは 1 solve 分 (本番 level 15 で MG 0.066 s、全 level 合わせて 0.5 s 程度)。`cg_auto_every=20` なら 1 step あたり 0.03 s。
  - 残る点: MG 選択中の level の CG の試行は MG の解からの warm start になる (その level が CG に切り替わった後の反復数とは少し違う)。

- **前処理 (5)**: level 15 の反復 26.9 → 11.1 (見込み 14)、1 反復は 0.58 → 0.97 ms (set-up を除くと 0.49 → 0.80 ms; ALLREDUCE の欄が
  0.22 → 0.37 ms に増えるのは遅いランクの掃引を待つ分)。level 15 の CG 10.8 → 8.3 s、poisson 10.3 → 8.6 s/step、step 時間 −3% で見積もりどおり。
  ただし `cg_auto_every=1` での比較。`cg_auto_every=20` と組み合わせた場合は未測定で、反復が同じ比 (0.41) で減るなら CG 18.5 → 14 s 前後
  (poisson 6.5 s/step、step 時間でさらに −2%)。
- poisson の CG 以外 (MG の level、levelmin、力) は 3 step で約 10 s (343085: poisson 28.8 − CG 18.5)。
- `sf - reduce` が 8.3 s/step (13%) で hydro - godunov に次ぐ。`star formation` (計算) は 3.6、`sf - ghostzones` は 0.3 s/step。2d-4 の予想どおり
  star formation の大半は集団通信の待ち (ランク間で 5–13 s/step とばらつく)。
- **本番の設定 (提案)**: 896x4、`KmpBlockTime "1"`、`omp_poisson=.true.`、`cg_auto=.true.`、`cg_auto_every=20`、`cg_warm_start=.true.`、
  `cg_precond` は使わない (得が 2–3% で実験的)。変更なしの 70.9 s/step に対して 62.7 s/step (−12%)。

### 2d-7. 残りの逐次部分のスレッド化: flag / rho / star formation (2026-11-20、クラスターの測定は末尾)

動機 (2d-4e の 343082 と 343087 の step 3): スレッド化済みの 4 行 (hydro・rt の godunov、cooling、poisson) は 896x4 と 512x7 で同じ時間、
通信・待ちは 512x7 が約 4 s 少ないが、未スレッド化の行が約 4.7 s 多くて引き分けになっている。待ちの原因は level ごとの grid 数の偏り
(合計の仕事量は max/mean 1.04 なのに level 15 は 1.38、level 14 は 1.75 で互いに逆相関 −0.91; level ごとの最大の和は平均の 1.83 倍) で、
単一の分割では消せない。ランクを減らす (1 ランクの領域を広げる) のが効くので、そのために逐次部分をスレッド化する。

- **粗い見積もりの訂正**: 896x4 → 512x7 の増え方から「grid 数に比例する計算」と「比例しない分 (通信・同期)」に分けると、896x4 での計算分は
  flag 0.8 / rho 1.0 / star formation 2.0 / refine 0.9 s/step 程度 (合計 4.7 s、step の 7%) で、残りは通信。スレッド化の得は 896x4 で 5% 前後、
  512x7 で 10% 前後の見込み (会話中に出した「896x4 で 57 s、512x7 で 50 s」は行の時間をすべて計算とみなした過大な値)。
- **スイッチ** (`&RUN_PARAMS`、既定 .true.、OPENMP=1 のビルドでだけ意味を持つ): `omp_flag`、`omp_rho`、`omp_sf`。どれも 1 ランクのその level の
  grid が `omp_mg_mingrid` (256) 以上のときだけスレッド化する。結果はスレッド数によらず逐次と bit 一致。
- **flag** (`omp_flag`): `hydro_flag` (上流の `hydro/hydro_flag.f90` を patch にコピー、`save` を外して grid の束のループを並列化) と
  `smooth_fine` の 3 つのループ (`count_nbors` の `save` も外した)。どの束も自分のセルの `flag1` / `flag2` だけを書く。`nflag` は整数の和。
  若い星まわりの refinement (`userflag_fine` の小野さんの追加部分) は隣の grid のセルに書くので逐次のまま。`init_flag` / `test_flag` /
  `poisson_refine` / `ensure_ref_rules` も逐次 (手元では flag の 1 割強)。`rt_hydro_flag` は本番で使っていない (`rt_refine` 未設定) ので未対応。
- **rho** (`omp_rho`、`rho_fine.f90` と `to_dice/rho_fine.f90` の両方; `multipole_fine` / `cic_from_multipole` / `cic_cell` / `rho_from_current_level` は
  両ファイルで同じ文面、`cic_amr` は DICE 固有の行だけ違う): `multipole_fine` は束ごとに自分のセルの `unew` だけを書くのでそのまま並列。
  隣のセルに足し込む CIC (ガスの `cic_from_multipole` と粒子の `rho_from_current_level`、2026-11-20 に粒子も対応) は次の形:
  - **スレッド化しないとき** (1–2 スレッド、`omp_rho=.false.`、grid が少ない level、OPENMP=0): `cic_cell` / `cic_amr` に `direct=.true.` を渡し、
    元のコードと同じくその場で `rho` / `rho_top` / `phi` に足す。
  - **3 スレッド以上**: 束を「スロット」(ガス 4 束、粒子 8 束) に分け、スロットを並列に計算する。`cic_cell` / `cic_amr` (`direct=.false.`) は足す先の
    セル番号と値を返すだけで、呼び出し側が元の加算順のまま「(足す先のセルの担当スレッド, スロット)」ごとのリストに書き足す。その後、各スレッドが
    自分の担当セル (セル番号 64 個ごとに同じ担当) の分だけをスロット順に足す。どのセルも逐次と同じ順序で項を受け取るので、スレッド数によらず bit 一致。
    levelmin の `multipole` は逐次。粒子の束を作る連結リストの走査も逐次。
  - 経緯: (1) 足し込みだけ逐次に残す 2 段 (計算 → バッファ → 逐次に足す) はガスで 2 倍止まり、粒子では 1 スレッドが 1.33 → 1.84 s と遅くなり 4 スレッドでも
    0.99 s (元の粒子の CIC は 1 粒子 100 ns 程度と軽く、バッファ経由の 1 往復で仕事が 5 割増える。`save` を外したせいではない)。
    (2) 「各スレッドがバッファ全体を走査して担当セルだけ足す」は逐次より遅い。(3) 担当別リストは、リスト長のカウンタを共有配列に置くと偽共有で
    2 スレッド 1.89 s / 4 スレッド 1.35 s; カウンタをスレッド私有にし、担当の計算から除算を除いて今の形。
  - 手元 (4 ランク、粒子が rho の 8 割を占めるテスト、rho の行全体): 変更前 1.33 s、1 スレッド 1.38、2 スレッド 1.39 (直接の経路)、4 スレッド 0.93。
    **4 スレッドで 1.4 倍止まり** (リスト経由の上乗せ 5 割と、逐次の走査・通信が残る)。本番はガスが主のはずなので比率は違う (未測定)。
  - 作業メモリ: リストは (スロットあたりの項数) × スレッド数² で、ガス 8192 項・粒子 2048 項 (phi は 2 倍) × 12 byte。7 スレッドで約 6 MB、14 スレッドで約 25 MB。
- **star formation** (`omp_sf`、`star_formation.f90`; `star_formation_fg.f90` は未対応): 保存量 ↔ 原始量の変換の 2 つのループを並列化。
  主ループは 2 段にした: 各セルの判定 `ok` と `sfr_ff` を先に全セル分求め (sf_model 6 / 7 かつ `sf_virial=.false.` のとき並列; それ以外は
  `getnbor` が `save` を持つので逐次)、Poisson の乱数はその後で元の順序のまま逐次に引く。乱数列が変わらないので結果は変更前と同じ。
  作業配列はセルあたり実数 1 + 論理 1。星粒子を作るループは逐次のまま。
- **refine は未着手**: 手元の内訳は ghost 交換 (全変数、level ごとに 5–6 本) 55%、`refine_fine` 41%、`build_comm` 4% で、本番でも
  ランク間でほぼ一様 (3.97–4.50 s/step) なので大半は通信と同期 (`refine_fine` は level ごとに `MPI_ALLREDUCE` 4 本)。スレッド化より前に内訳が要るので
  `amr_step` に timer ラベル **`refine - comm`** (`build_comm` + 交換) を足した。**`refine` の行は今後 `refine_fine` だけ**になる
  (以前の `refine` = 新しい `refine` + `refine - comm`)。
- 検証 (`run/a_warm4` の namelist を 10 step、4 ランク、level 6–8; 変更前 = HEAD 30db6e70 の worktree): OPENMP=1 + BITREPRO=1 の 1 / 4 スレッドは
  変更前の 1 スレッドと 80 ファイル bit 一致 (`th_b1` / `th_r1` / `th_r4`)。星ができる設定 (sf_model=1、64 個) でも 1 / 3 スレッドで 60 ファイル一致
  (`th_m1b` / `th_m1n1` / `th_m1n4`)。OPENMP=0 のビルドも変更前と 80 ファイル一致 (`rs_D` / `th_noomp`)。sf_model=6 はこの初期条件では
  条件を通るセルの `sfr_ff` がすべて 0 で星ができないので、デバッグ用ビルドで `sigs + scrit` を混ぜたチェックサムを取り、1 / 4 スレッドで
  45704 セル分が bit 一致することを確かめた (**sf_model=6 で実際に星ができるランでの確認はしていない**)。
- ks 設定 (OPENMP=1、DICE なし = `rho_fine.f90` の方、`sf_fg`; `run/kstest_ts` の restart を 12 ランク 31 step、`th_ksb` / `th_ks1` / `th_ks2`):
  1 スレッドは変更前のバイナリと 72 ファイル bit 一致。2 スレッドと 3 スレッド (粒子の CIC のリスト経路を通る、`th_ks3`) は
  `omp_hydro/omp_rt/omp_cooling/omp_poisson=.false.` (新しい 3 つだけスレッド化) で 72 ファイル一致。全部スレッド化した 2 スレッドは 60 ファイル不一致だが、これは既存の成分がスレッド数で丸めを変えるため (今回の変更ではない)。
  この zoom ランでは 2 スレッドでも flag 9.3 → 8.8 s、rho 18.5 → 19.1 s と縮まない (1 ランクの level あたりの grid が少なく、rho は粒子が主)。
- 手元の時間 (4 ランク、他の成分もスレッド化、1 → 4 スレッド): `hydro_flag` 0.345 → 0.121 s、`smooth_fine` 0.192 → 0.082、flag 全体 0.44 → 0.16;
  `multipole_fine` 約 2 倍、rho 全体 1.38 → 0.93; star formation 0.110 → 0.047 (同じ機械で hydro - godunov は 2.8 倍)。
  他の成分を逐次にして測るとスレッドが眠っていて伸びが小さく見える (flag 0.43 → 0.32) ので、測るときは全部スレッド化した状態で。
- **クラスターで見るもの**: 896x4 と 512x7 (どちらも warm + `cg_auto_every=20`) で `flag`、`rho`、`star formation`、`refine`、`refine - comm` の行と step 時間。
  512x7 が勝つなら 256x14 も試す価値がある。`refine` (= `refine_fine`) が大きければ `authorize_fine` と走査のループが次の対象。


#### 2d-7 のクラスター測定 (2026-11-20 夜、バイナリ c4ce85ac)

ログは 2d-4e と同じ場所。3 本とも level 15 の RT subcycle は step 2 / 3 が 256 / 256。数字は s/step (step 2–3 のランク平均)。

| ジョブ | 構成・設定 | step 2 / 3 (s) | 比較相手 | poisson | flag | rho | star formation | refine (+ comm) | sf - reduce |
|---|---|---|---|---|---|---|---|---|---|
| 344393 | 896x4、warm + `cg_auto_every=20` | 61.5 / 60.8 | 343085: 63.8 / 61.5 | 8.3 (7.9) | 1.25 (1.85) | 2.30 (2.24) | 2.07 (3.60) | 1.65 + 2.63 (4.3) | 8.3 (8.4) |
| 344397 | 896x4、warm + `cg_auto_every=1` | 62.0 / 59.8 | 343082: 65.2 / 64.8 | 8.8 (10.3) | 1.26 (1.82) | 2.27 (2.27) | 2.04 (3.61) | 1.66 + 2.62 (4.2) | 7.9 (8.3) |
| 344396 | 512x7、warm + `cg_auto_every=20` + `omp_cg` | 57.8 / 58.1 | 343087 の step 3: 65.4 (`cg_auto_every=1`) | 7.6 (10.4) | 1.35 (2.42) | 2.55 (2.98) | 2.07 (5.13) | 1.98 + 2.98 (5.0) | 6.3 (7.0) |

括弧内は比較相手の値。

- **512x7 が 896x4 に勝った: 57.9 対 61.2 s/step (−5%)。** 前回は引き分け (65.4 対 64.8)。本番は 512x7 に変える価値がある。
- スレッド化の得 (896x4): star formation −1.5、flag −0.6、rho は変わらず (+0.05)、合計 −2.1 s/step (step 時間は 62.7 → 61.2、−2.4%; 見込みの 5% の半分)。
  512x7: star formation −3.1、flag −1.1、rho −0.4 で −4.6 s/step。343087 との差 7.3 s の残りは `cg_auto_every` 1 → 20 の poisson (−2.8)。
- **rho は効いていない** (896x4 で 0、512x7 で −14%)。行の大半が通信か、粒子の CIC のリスト経由の上乗せで相殺されている。内訳は未調査。
- **refine の計算 (`refine_fine`) は 1.65 s/step (512x7 で 2.0)、通信が 2.6 (3.0)。** スレッド化しても得は最大 1.5 s 程度。
- **`cg_auto` の解き直し (30db6e70) は効いた**: `cg_auto_every=1` の poisson が 10.3 → 8.8 s/step になり、`cg_auto_every=20` (8.3) との差は試行の時間だけになった。
  step 時間は 344397 と 344393 で同じ (60.9 と 61.2)。level 7–10 は 3 本とも MG のまま。
- `feedback` が 1.29 → 1.49 (896x4)、0.76 → 0.87 (512x7) と増えた。原因は未調査 (前の行の待ちが移っただけの可能性)。
- 残っている大きい行 (512x7): hydro - godunov 9.2、poisson 7.6、sf - reduce 6.3、hydro - rev ghostzones 5.3、cooling 3.5、refine - comm 3.0、rt - godunov 3.0。

#### 256x14 (344731、2026-11-22、`hu_lv15_omp14/`、namelist は 344396 と同じ)

32 ノード × 8 task × 14 スレッド (NUMA node 1 つに 1 task)。`output_00051` は Mera で 256 ランクに書き直したもの。

- 壁時計 85.7 / 84.3 / 74.1 s。ただし **level 15 の RT subcycle が 486 / 512 / 394 回** (344396 は 511 / 256 / 256) で、そのままは比べられない。
- RT に比例する行 (rt - godunov、cooling、rt - ghostzones、rt - rev ghostzones、rt - set unew、rt - upload fine) は step 2 / 3 で
  25.6 / 19.2 s、subcycle あたり 0.050 / 0.049 s (512x7 は 11.2 s、0.044 s)。256 回に揃えると 12.5 s。
- RT 以外は 58.3 / 54.5 s (512x7 は 46.5 s)。**RT を揃えた 1 step は 67–71 s 相当で、512x7 (57.9 s) より 16–22% 遅い。**
- 遅くなった行 (512x7 → 256x14、s/step): poisson 7.6 → 11.5 (level 15 の CG 反復が平均 17.8 → 22.0 回、1 解法 0.022 → 0.029 s)、
  refine 2.0 → 2.9、refine - comm 3.0 → 4.2、rho 2.5 → 3.6、courant 1.45 → 2.3、hydro - set uold 0.9 → 1.9、hydro upload fine 0.35 → 0.67、
  flag 1.35 → 1.65、star formation 2.07 → 2.40。逐次のまま残っている部分が 1 ランクの grid 数 (2 倍) に比例して増え、スレッド化した部分は
  14 スレッドで伸びきらない。速くなったのは feedback (0.87 → 0.50) と hydro - rev ghostzones (5.3 → 4.9) だけ。hydro - godunov は 9.2 → 9.5。
- **結論: 256x14 は採らない。本番は 512x7。**

### 2d-5. 512 ランク × 7 スレッド (flat-339876、2026-11-15、変更前のバイナリ、`KMP_BLOCKTIME=0`)

`hu_lv15_test_512/output_00051` (Mera で 512 ランクに書き直し) から 3 step、32 ノード × 16 task × 7 スレッド。比較は 896 × 4 の 339731。
- 配置は意図どおり: task 0 がコア 0-6、task 2 がコア 7-13 (14 コアの NUMA node に 2 task、task は socket 交互)。
  ノード内の位置別の hydro-godunov は 16 位置とも 10.0-10.1 s で、339731 の位置 6/7/20/21 の遅れ (10.8/11.4/10.8/10.7 s) は消えた。
  ランク間の max/min も 1.50 → 1.31。
- coarse step の壁時計は 97.3 / 86.0 / 85.8 s (339731 は 104.0 / 90.9 / 88.8 s) だが、**RT の履歴が違うので合計では比べられない**:
  339876 は dt<1e-5 の fine step が多く (193 / 236 / 43; 339731 は 0)、level 15 の RT subcycle が 1 step 256-268 回 (339731 は 512 回)。
- step 2 のランク平均 (339731 → 339876): poisson 18.9 → 28.0 s (逐次、1 ランクの grid が 1.75 倍; CG level 15 は 1 解法 0.025 → 0.040 s)、
  hydro-godunov 9.8 → 10.1 (コアあたりの効率は同じ)、rt-godunov 6.1 → 3.3、cooling 7.1 → 3.7 (subcycle あたりでは +8% / +4%)、
  star formation 15.5 → 12.7、hydro-rev ghostzones 6.8 → 5.5、rt-ghostzones 6.8 → 2.9 (subcycle あたり -15%)、
  refine 4.6 → 5.1、rho 2.2 → 3.0、flag 1.7 → 2.6、courant 1.1 → 1.4。
- RT の仕事を 512 subcycle に揃えると 339876 は 1 step 約 98 s 相当で、**変更前のバイナリでは 896 × 4 のほうが 7% ほど速い**。
  差は逐次の poisson (+9 s) と grid 数に比例する逐次部分 (refine/rho/flag/courant +2.4 s) で、待ち (star formation、ghostzones) は減っている。
- poisson をスレッド化した新バイナリ (2d-4) では逆転しうる: 手元の伸び (4 スレッド 2.4 倍、7 スレッド 2.8 倍) がそのまま出れば
  どちらも 1 step 82 s 前後。**4 と 7 のどちらが速いかは新バイナリで両方の配置を測って決める** (判定なし、`KMP_BLOCKTIME=1`)。
- cg_auto の選択は 18 CG / 8 MG (339731 は 16 / 10) で、level 10 以上が CG、9 以下が MG という分かれ方は同じ。

### 2d-6. SN による dt の急落 (2026-11-16 に気づいた、未対応)

本番 (`hu_lv15`) では、ガスがほとんど無いセル (0.002-0.07 Msun) に SN の熱が入って `FEEDBACK T2` = 5e10-2.5e12 K になり、dt が
1/30-1/70 に急落してからゆっくり戻る。3 step のテストでは毎回 1-2 回起き、coarse step 1 つで進む時間が半分以下になる
(339731 は 0 回で 0.31 Myr、339876 は 2 回で 0.14 Myr)。RT の subcycle 数もこれで変わる (dt が小さいと 1 fine step に 1 回)。
SN は乱数で起きるのでランク数・スレッド数・丸めで履歴が変わり、**3 step の合計時間はランごとに比べられない** (項目別・step 別に見る)。
- 既存ログの粗い集計 (`hu_lv15/flat-*.out` の最後の 20 本、114 Myr、47 万 fine step): 急落 (前の step の 0.3 倍未満) は 158 回
  (coarse step 11 個に 1 回)。急落前の dt の 0.8 倍に戻るまでを「回復中」と数えると fine step の 83% がそれに当たり、急落が無ければ
  同じ時間を 4 割の step 数で進める計算になる。回復の定義に強く依存する粗い数字 (急落前の dt がたまたま大きかった場合は過大) だが、
  スレッド化より大きい損失の可能性がある。
- ks の SN reservoir (冷却半径を分解できる球に熱を入れる; `SN_reservoir`、移植済み) は、まさにこの dt の急落と overcooling の両方を防ぐために
  作られた実装 (杉村さん → ユーザーの説明、2026-11-16)。対処するならこれを使う (温度上限で隣に配る等の簡略版を新しく書く理由は無い)。
  ただし SN の入れ方が変わるので前回のランとの比較には向かない (ユーザー判断 2026-11-16: 300 Myr からの流し直しでは変えない)。
  順序: 速い OpenMP の設定を決める → 物理を変えずに流し直し → そのあと `SN_reservoir=.true.` を本番 snapshot から 3 step 試して dt の履歴を見る。調べるなら、急落を起こした SN のセルの密度・温度の分布と、
  dt を決めているセルの数 (1 セルか bubble 全体か) から。

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

## 3b. patch/ks からの移植 (merge_ks ブランチ、docs/merge_ks.md)

- `test_flag` の修正 (子の refine を所有者から配る) は ks と同じく常時有効にした (2026-11-13; 単一分割でも
  境界例で孤児 grid を作るバグなので)。既定テストは bit 一致。本番では境界例で buffer oct が少し増える方向。
- Mera.jl の two-tier 出力 (`hilbert2t`) と QUADHILBERT 対応は後回し (merge_ks.md 末尾)。
- 星形成は移植済み (2026-11-13, `sf_fg`, `star_formation_fg.f90`; merge_ks.md「星形成の移植」)。残り:
  - ugvtk 出力は移植済み (2026-11-14; merge_ks.md「ugvtk 出力の移植」、主スイッチ `ugvtk`)。未移植: `fout_*` 系の付録出力。
    ks と同じ要確認点: `output_around_SFC_at_given_age` 経由の uold 往復変換 (丸め差)、`set_ufront(lev,3)` の無条件上書き。
    `my_utils` の `get_cpu_of_pos` は boxlen≠1 で grid 幅の単位が違っていたのを修正 (2026-11-14; ks は boxlen=1 で無影響)。
    年齢出力は 1 step に 1 つ (fine step より短い 10/100 kyr は飛ぶ、ks と同じ)。
  - `SFC_flg_thin` / `flg_thin_star_cell` の cooling 側、`fJeans_cool`。
  - Pop III (family 3) と SFC は takashi の feedback/RT では `is_star` でないので放射も SN も出ない → フィードバック段階。
  - ks の part ファイルからの restart: `mp0` の有無 (ks は star_by_star/fire3_yields のときだけ) と `msp0` (takashi のみ)。
  - FG の stdout (星 1 個ごとの print) の抑制。`random_value_0to1` の乱数状態は restart で保存されない (ks と同じ)。
  - `sf_fg` では `f_w>0` を禁止している (debris family なし)。sink は移植済み (下)。
  - 非局所枯渇の一括版 `deplete_gas_from_spheres_batch` (`nonlocal_depletion_batch`, 既定 .true.) を入れた
    (2026-11-13; merge_ks.md「非局所枯渇の一括処理」)。`nH_postSF>0` の `deplete_gas_with_fixed_postSF_density` と
    Pop III の `deplete_gas_around_pos` は ks の逐次版のまま (個数が少ない)。SN reservoir も同じ octree 走査
    (`sphere_walk`) が使える。
- フィードバックは移植済み (2026-11-13; merge_ks.md「フィードバックの移植」): Pop III SN、SFC→PSC、SN reservoir
  (`SN_reservoir`)、冷却停止タイマー (`fcool_stop`)、Pop III の放射。残り:
  - 未移植: mechanical/spread/nonlocal_SN、star-by-star の SN (Inoue+2011, FAM_DEAD)、`pop2_SN_model`、hypernova X 線 (`isXray`)。
  - reservoir の pending エネルギーを含む checkpoint からの restart は経路のみ確認 (pending が 0 の output からの restart は一致)。
  - 一括非局所枯渇と SN reservoir は別の octree 走査 (`sphere_walk` と `get_sn_sphere_stats_events`)。統合の余地。
- sink/BH は移植済み (2026-11-13; merge_ks.md「sink/BH の移植」、`sink_particle_fg.f90`)。残り:
  - 未移植: clump ベースの sink 生成 (`make_sink_from_clump`; ks 本番も bhpop3 のみ)、sub-grid AGN の種入れ
    (`sgagn_seed='mhalo_based'`, halo finder が要る)、旧 `accrete_sink` (AGN blast/`sink_drag`)、`check_particle_consistency` (ks でも呼び出しはコメントアウト)。yt が読む `sink_file_descriptor.txt` は出す。
  - sink の csv は上流の `output_sink_csv` (msmbh/dmfsink は 0; ks の 13 列は年齢が tsink になる以外すべて含む)。
    ks 形式の binary `sink_XXXXX.out<rank>` (全ランク同じ内容; yt は全ドメインのファイルの存在を要求するので rank 1 だけにはできない)
    と yt 用 descriptor、ks のテキスト表 `sink_XXXXX.info` (`output_sink_fg`) も出す。yt (dev 4.5 / 4.3.1) で `sink` と
    `sink_csv` の両方が読めることは確認済み (2026-11-14)。
  - `check_sink_vars_among_cpus_fg` は `sink_check_vars=.true.` のときだけ (sink 1 個・呼び出し 1 回あたり ALLREDUCE 7 本)。
  - `get_BHL_accretion_rate_fg` は非 cosmo でも G=1 で動くようにした (ks は停止)。`d_sink` 未指定でも止めない。
  - DICE の `star_as_nbody` 粒子は family 5 = `FAM_DEAD` になるので `sink_mass_min_pop2BH>=0` にすると全部 sink 化される。
- ks との比較ラン (2026-11-14; merge_ks.md「ks との比較ラン」): MUSIC zoom (haloD_l5_14_z9) 6 step で ks 版と初期状態が全桁一致、
  以後は密度 1e-6・運動量 3e-5 の相対差、grid 数一致。合わせるのに `omega_b` (ks は omega_b0 を無条件代入) と
  初期電離度 (Peebles の残留電離を `init_flow_fine` に移植、`rt_is_init_xion=.false.`) が必要だった。残り:
  - 長いラン (12 ランク、z=21 まで): 最初の Pop III 形成 (pop3_model=6、8 連星 2413 Msun)、3 つの直接崩壊、BH sink 3 個の
    生成と合体 (2120 Msun)、Bondi 降着まで ks と同じイベント・同じ数値で一致 (merge_ks.md「長いラン」)。その際 takashi の
    FIRE-3 feedback の年齢換算 (cosmo で a² 倍小さい) を修正。
  - PISN 以後 z=20 まで (2026-11-15; merge_ks.md「PISN 以後の比較」): Pop III の寿命末 7 件と SN 後の dt、メタル質量 (5e-5) は一致。
    HII 領域が takashi で 1-3 割小さいのは ks の試験 namelist の群 5 (既定値 13.6 eV〜上限なしの重複した電離群) が原因で、
    takashi 側に同じ群を入れると電離質量の差は 13% → 2.7%。杉村さんに伝える (`run_l5.sh` の `groupL0/L1` が 4 群ぶんしかない)。
    残りの 2.7% は `attn_after_chem` の解き方 (ks: 更新後の状態で吸収・収束判定なし) で、`attn_new_abund=.true.` + `rt_skip_convergence=.true.`
    にすると 0.4% まで合う。PISN の注入は一致、その後の高温ガスは takashi が速く冷える (冷却テーブルの差と推定、冷却率は未比較)。
    残り: 高温・高金属量ガスの冷却率の比較、Pop II 形成の比較 (z=20 では未形成)。
  - 収束判定 (2026-11-15; merge_ks.md「吸収の順序と収束判定」、`tests/ifront_thick.nml`): 厚いセルの R 型波面では `rt_skip_convergence=.true.`
    が通常の刻みで 4 割以上多く電離する (刻みを細かくすると判定ありの値に寄る)。DICE 本番 (`xd_rt.nml`) は `.true.`。本番の密度・解像度で
    どれだけ効くかは未評価。`namelist/ks/*.nml` は既定 (判定あり) のまま。コストは手元の DICE (fb6 相当 24 step、4 ランク) で
    cooling 105.5 → 167.3 s (+59%)、全体 +26%。このテストはガスがほぼ全部電離していて (HI 2e3 / 1.3e8 Msun) 電離波面が無く、結果への
    影響は測れなかった (HeIII の質量差 0.8%)。本番 snapshot から `rt_skip_convergence=.false.` で数 step 流して電離質量と時間を比べる。
    本番の出力での見積もり (`tests/tau_step_stats.py`、`hu_lv15/output_00051`、3650 万セル): 電離光子の吸収のうち RT 1 step の光学的厚さが
    1 を超えるセルの分は 9.4% (3 超は 2.1%)、0.3 未満が 72%。吸収の 63% は光子が豊富なセル (a>=1、再結合との釣り合いで誤差なし)、
    光子の乏しい波面側 (a<1) は 37% で吸収加重の平均 tau_step は 0.4-0.6。判定なしの電離しすぎは全体で高々 1 割程度と見積もれる
    (電離波面テストの 4 割増の条件 = tau_step~4 に当たる吸収は 2%)。直接の確認は本番の 3 step テスト (nstepmax=7578, foutput=7578)。
    ただしこの見積もりは判定なしで作った snapshot の状態から数えているので偏る: 厚いセルは誤差で早く電離されてしまい、統計に残らない。
    密度の高いガスに埋もれた光源 (解像されない Stromgren 球; `ifront_thick.nml` を nH=30、0.1 Myr = 再結合時間の 25 倍に変えたもの) では、
    電離原子数は判定あり 1.397e60 (ks の順序も同じ、上流の既定 1.368e60、刻み 1/8 で 1.411e60) に対し判定なし 1.678e60 (+20%)、
    飛行中の光子は 12 倍。平衡でも判定なしは電離しすぎ、閉じ込められるはずの光子が漏れる方向。本番の 3 step テストでは
    nH>1 のガスの電離質量が判定ありでどれだけ減るかを見る (再結合時間が短いので 3 step で落ち着く)。
    **本番の 3 step テスト (2026-11-16、`hu_lv15_conv` 340010 / `hu_lv15_noconv` 340011、896 × 4、`new_dec` f5ba86de のバイナリ)**:
    2 本は SN の履歴が同じ (同じ fine step で dt 急落、終了時刻 t=33.562912 = 0.18 Myr 後) で直接比べられる。
    コスト: coarse step 2/3 が 74.3/73.2 → 77.7/77.2 s (+5%)、cooling 3.5 → 3.9 s/step、Coolstats の平均ループ 1.004 → 1.206 (最大 681 → 3289)、
    増分の半分以上は star formation の行 (待ち、12.0 → 14.5 s)。結果 (`tests/ion_by_density.py`、判定あり / なし - 1):
    HII は nH 1-10 で -0.5%、10-100 で -5.8%、100 以上で -6.4% だが、そこにある HII は 3.7e4 Msun (その密度の水素の 0.25%) しかなく、
    HI は +0.004〜+0.07%、H2 は ±0.02% 以内、全体の HII は -0.005%、熱エネルギーは nH>10 で -2〜3%、飛行中の電離光子は -0.07〜-0.14%。
    **判定なしでも本番の結果はほとんど変わらない** (電離しすぎは密度の高いガスの中の電離成分 6% にとどまり、漏れ出す光子も 0.1%)。
    `xd_rt.nml` は `rt_skip_convergence=.true.` のままでよい。
  - FIRE-3 の年齢の評価点を ks 式 (step 中央、最初の step は年齢 0 からの積分) に揃えた (2026-11-14)。DICE の結果は
    変わる (SN 率の評価点が半 step ずれる) ので、以後の回帰の基準はこの版の出力。
  - 初期温度 0.5% (`rt_evol_single_cell` の冷却率・定数)、`x_MIN` (ks 1e-12 / takashi 1e-20) の差は残したまま。
  - 未移植で ks 側を off にして比べたもの: `flg_voffset_ini` (初速オフセット)、`fof_find`、`flg_thin_star_cell`、`pop2_SN_model`。

## 3c. sf_model=6 の σ の定義 (2026-11-24、docs/sf_model6_sigma.md)

- 本番の output_00038 で、星ができるセルの σ は主に収縮流を測っていた (星形成率の 99.8% が `div v < 0`)。
- 足したもの: `sf_sigma_noconv`、`sf_alpha_meandens` (既定 .false.、既定では bit 一致)。`sf_b_turb` の既定値を 0.4 → 0.4√3 に変更
  (**namelist に書いていないランは結果が変わる**; 以前の値は `sf_b_turb=0.40000000596046448`)。
- やり残し:
  - 本番の分解能で流していない。output 38 から 20-30 Myr 流して、星形成率・星のできる密度・SN の dt への影響を見る。
    瞬間の見積もりでは 2 つのスイッチ + 新しい b で星形成率 28 倍 (自己調整でどこまで戻るかは不明)。
  - father セルが入る箱 (level 境界) での `sf_sigma_noconv` は数値を確かめていない。
  - 乱流エネルギー K の方程式 (KT20 式 4-7) は未実装。
  - `sf_virial=.true.` の経路と `sf_model=8` には入れていない。

## 3d. ダスト (放射) 温度 TR と IR 群の判定 (2026-11-24、docs/dust_temperature.md)

- `is_kIR_T` の TR を釣り合いの式で解くようにした (上流の 1 回の式は UV 支配の所で TR が 2-4 倍高い)。`rt_isIR_alt` の
  IR 群の subcycle 判定は再放出を足した正味の変化で行う。旧来の式は残していない。DICE の本番では結果は実質変わらない。
  IR 群へ移す量 (UV・可視のダスト吸収分と IR の再放出) は吸収後の光子数で数える (以前はガスが吸収した光子もダストの分に数えていた)。
- やり残し: 群 1 に星の近赤外が入る、TR に CMB の下限が無い。本番の分解能・`rt_skip_convergence=.false.` での計算時間は未測定。

## 4. 運用上のメモ

- **restart ファイルの保護 (2026-11-16、patch/ks から移植)**: `dump_restart` (`output_amr.f90`) は `output_restart.tmp/` に書き、全ランクが
  書き終えてから `output_restart/` に rename する。前の組は次の dump が始まるまで `output_restart.prev/` に残る。dump の途中で job が
  落ちても `output_restart/` は常に完全な組。ディレクトリ操作は `posix_fs.f90` (libc を ISO_C_BINDING で呼ぶ; patch/ks からコピー、
  `bin/Makefile*` 5 本の MODOBJ に `posix_fs.o` を追加)。読む側 (`nrestart=-1` → `output_restart/`) は変えていない。
  確認 (dice_rt_amr、4 ランク、`restart_hrs=0.001`): dump 4 回で `output_restart/` (nstep 7) と `.prev/` (nstep 5) が残り `.tmp/` は無い、
  `nrestart=-1` で再開して完走。再開後の出力は通しのランと bit では一致しない (16/16 ファイル) が、差は丸めの大きさ
  (総和・min/max の相対差が最大 1e-9、大半は 1e-13 以下、セル数と log のエネルギーは同じ): restart ファイルは倍精度で値をそのまま持つので、
  再開後に演算の順序が変わる (粒子リストの並び、再開直後の load balance など) ためと考えられる。ディスクには書き込み中に最大 2 組 (本番で 24 GB × 2) が載る。

- **2 段の出力ディレクトリ `output_XXXXX/group_YYYYY/` (2026-11-22 に修復)**: `IOGROUPSIZEREP` は namelist の変数ではなく
  `amr_parameters.f90` の値 (既定 0 = 1 段、本番はこれ)。N にすると N ランクごとに `group_YYYYY/` に分ける。上流では 2019 年の
  9988a054 (ディレクトリ作成を cpu 1 だけに) 以来、cpu 1 のグループしか作られず壊れていた。直した点: (1) `create_output_dirs` は
  cpu 1 が上の段を作り、2 段のときは全ランクが自分のグループを作る (NOSYSTEM の `filedir(1:13)` 決め打ちもやめた: 空の
  `output_restar/` ができていた原因)。(2) `read_params` (両方のコピー) の restart 用 info の存在確認が上の段しか見ていなかったので
  `group_00001/` を見る。確認 (gfortran、`IOGROUPSIZEREP=2` にしたビルド、DICE RT 4 ランク): 書き出しが 1 段のランと bit 一致 (24 ファイル)、
  `output_restart/group_*/` から `nrestart=-1` で再開した結果が、同じ dump を 1 段に並べ直して再開した結果と bit 一致 (44 ファイル)、
  `nrestart=2` も同様 (24 ファイル)。既定 (1 段) は変更前と bit 一致。**未確認**: Intel (NOSYSTEM、`PXFMKDIR`) 側、sink と ugvtk を
  使うランの 2 段出力、Mera.jl / yt が 2 段の output を読めるか。

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
- 未追跡ファイルは 2026-11-15 に整理した: `*~`、findent 整形前の `backup_*` (git 履歴と空白以外同一)、
  `rt_cooling_module.f90.back`、`to_dice/init_flow_fine_corrected.f90` (a6263c68 の `init_flow_fine.f90` と同一)、
  `to_dice/format_fortran_script*.sh` (`patch/takashi/` の追跡済みと同一)、`*.nml.bak` を削除。
  `.gitignore` に `*~`, `*.bak`, `patch/takashi/logs/` (本番ログ、未追跡のまま置く), `.claude/settings.local.json` を追加。
- 2026-10-29: per-level 分割のために patch に置いていた `load_balance.f90, refine_utils.f90,
  sink_particle.f90, flag_formation_sites.f90, init_tracer.f90` は削除した (上流が使われる)。
