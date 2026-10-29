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

### 1d. balance_per_level で本番 896 ランクの grid 所有者と authorize が食い違う — **未解決**
- 2026-10-27: `make_grid_fine` が新 grid のセルの `cpu_map` を `bound_key(:,ilevel)` で決めていた
  (正しくは `bound_key(:,ilevel+1)`) のを修正。`build_comm` に「他 CPU が自分の仮想 grid として
  持っているのに自分に実体がない」検出を追加。
- 2026-10-28: 修正版でも本番 (flat-329902, `balance_level_min=13`) で 5 coarse step 目に level 15 で
  検出が発火 (cpu 622, 27 grids)。手元 (8/32 ランク, levelmax 8) では再現せず、コード上の矛盾も
  見つかっていない。`build_comm` の検出を集団診断 (`debug_missing_owner`: holder 側の親セル情報と
  owner 側の木の走査を出力) に拡張した。次の本番ジョブのログ (`MISSING` / `OWNER` 行) で原因を
  特定する。詳細は `docs/level_domain_decomposition.md` 4b-2。
- 疑うべき候補: 所有者側が親セルを持っていない (祖先のどこかで authorize が外れている) か、
  持っているのに refine していない (flag1 の同期漏れ) か。`OWNER` 行の最深 level と T/F で切り分ける。
- 2026-10-28 (flat-330537, 誤って per-level で走ったジョブ): 診断が初めて出た。cpu 660 / level 15 /
  27 grids (holder 664, 665)。owner 660 は level 13 まで持ち level 14 の grid が無い (`son=0`)。
  holder 側の親セルは `cpu_map=0, son=0`、親 grid 所有者 5242 (> ncpu) と不正で、stale な仮想 grid
  か `father` の破損を疑う。docs 4b-2 末尾に詳細と再開時の次の一手。棚上げのまま。

### 1e. (修正済 2026-10-28) `feedback.f90` の SN イベントログで添字違い (上流由来)
- `feedbk` の `sf_log_properties` ブロック末尾 `typep(ind_part(i))%tag` が、ループ変数 `j` の
  ところで古い `i` を使っていた (上流 `pm/feedback.f90:405` も同じ)。`sf_log_properties=.true.`
  で SN が起きると `typep(0)` を読み、-O2 では segfault、bounds check では
  "Index '0' of dimension 1 of array 'typep'"。`j` に修正。`patch/mizuki/feedback.f90` にも同じ行がある
  (そちらは触っていない)。本番 namelist は `sf_log_properties=.true.` なので本番でも踏み得た。
- 2026-10-29: origin/mizuki acb88172 でログ出力ブロックが MHD のループ注入の前に移動した (こちらにも
  移植済み、hydro ビルドでは bit 一致)。mizuki 側はまだ `typep(ind_part(i))` のまま。小野さんに連絡要。

## 2. level ごとの domain decomposition (balance_per_level) の改善案

設計と検証結果は `docs/level_domain_decomposition.md`。

**2026-10-28 の結論: 見送り。** 本番 (896 ランク) の 4 ジョブで、RT を除いた step 時間は単一分割と
差がなく (165-177 s vs 170-188 s)、level 15 を完全均等にしても RT subcycle 1 回 0.39 s は変わらず、
grid メモリは 27% → 36-39% に増えた (docs 4c)。加えて原因未特定の不整合 (1d) が残る。
既定 `balance_per_level=.false.` なら単一分割と bit 一致なのでコードは残す。手元の長い再現テスト
(levelmax 10, 24 ランク, level_min=8, τ=0) は 121 coarse step まで不整合なしで、途中で止めた。
次に効くのは RT step の内訳 (通信か冷却か): 単一分割で `nstepmax` を小さくして正常終了させ timer 表を見る。

- **ghost 抑制 (2026-10-27 実装済)**: `balance_level_min` (粗い level は指定 level の分割をコピー) と
  `balance_tolerance` (許容不均衡の範囲で境界を 1 つ細かい level に寄せる)。設計と目安は
  `docs/level_domain_decomposition.md` 4c。**本番 (896 ランク) での値探索は未着手**:
  単一分割 / per-level τ=0 / τ=0.3 / τ=0.45 / balance_level_min=13 を同じ restart から流し、
  ghost 合計・`Level l has N grids (min,max,avg)`・coarse step 時間を比べる。
  これでも ghost が 2 倍以内にならなければ、項 1 (自分の grid の子が他人の grid) を消すには
  upload を子の所有者側で行う大改造が要るので、per-level 自体を見送る判断もあり得る。
- 4 ランク DICE+RT テスト (`tests/dice_rt_amr.nml`) の timer では、load balance 0.97→0.88 s,
  godunov の CPU 間ばらつき 0.015→0.031 s, 合計 9.39→9.23 s で差は小さい (4 ランクでは
  偏り自体が小さいため)。効果の評価は本番の CPU 数でないと意味がない。
- **本番規模での実測**: `Level l has N ghost octs (min, max, ghost/active)` の行を
  level ごとに見て、ghost とメモリ (`mem=`) と `timer` の通信項目を単一分割と比較する。
  4 ランクのテストでは領域のずれが極端に出るので、本番の CPU 数での評価が必要。
- oct 数が ncpu 以下の level では 1 CPU あたり 0-1 oct になる。問題になるなら
  そういう level は粗い level の分割をコピーする分岐を足す。
- 後処理ツール: level ごと分割時は `info_XXXXX.txt` の DOMAIN 表に全域を書くため、
  Hilbert key でファイルを選ぶツール (yt, pynbody, utils/f90/*) は全ファイルを読む
  (遅いが正しい)。必要なら `info_levelkeys_XXXXX.txt` を読むように改造する。
- 旧コードは新しい `amr_*.out` (末尾に level ごとの key が追加) を restart できない
  (level 1 の分割だけ読まれ cpu_map と矛盾する)。

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
  `balance_per_level=.true.` でも流せば (abort は 4 step 目) per-level の表も直接比較できる。
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

## 3. 通信の隠蔽 (Phase 2) のやり残し

- 実装済み: 全変数を 1 メッセージに集約 (`make_virtual_fine_dp_nv`)、`unew`/`rtunew`
  の reverse 通信を godunov と重ねる非同期化 (`make_virtual_reverse_dp_start/end`)。
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

- `patch/takashi/` には未追跡ファイルが残っている (`flag_utils.f90`,
  `rt_cooling_module.back.f90`, `to_dice/init_flow_fine_corrected.f90`,
  `backup_*` ディレクトリ, `format_fortran_script*.sh`)。必要なものは `git add`、
  不要なものは削除して状態をはっきりさせたい。
- 今回追加したファイル (要 `git add`): `load_balance.f90, refine_utils.f90,
  sink_particle.f90, flag_formation_sites.f90, init_tracer.f90,
  docs/level_domain_decomposition.md, TODO.md, tests/` (2026-10-21 commit 済み)。
