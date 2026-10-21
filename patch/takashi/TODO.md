# やり残し・既知の問題 (引き継ぎメモ)

最終更新: 2026-10-21 (level ごとの domain decomposition の作業時点)

## 1. 既知のバグ

### 1a. `patch/takashi/flag_utils.f90.mizuki` (git 未追跡, ビルド対象外に改名済み) の若い星まわり refinement — **要修正**
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

### 1b. `to_dice/init_part.f90`: `zp` の確保条件
- `zp` (粒子金属量) は `if(star.or.sink)` の中で `if(metal)` のときだけ確保されるが、
  `to_dice/init_flow_fine.f90` は `metal` だけ見て `zp` を読む。`eps_star=0` で
  `metal=.true.` だと未確保配列を読む (本番 namelist では起きない)。
  `if(metal)` の確保を `star.or.sink` の外に出せばよい。

### 1c. `rt_init.f90` ~316 行: `rt=.false.` かつ `neq_chem=.false.` で落ちる
- NIONS=4 ビルドで ion の添字 (`ixHII` 等) が設定されず `group_csn(…,0)` を参照。
  本番では常に rt=.true. なので実害なし。テスト時の罠。

## 2. level ごとの domain decomposition (balance_per_level) の改善案

設計と検証結果は `docs/level_domain_decomposition.md`。

- **ghost oct の増加を抑える重み付け** (優先度高):
  level `l` の分割を決めるコストに、細かい level `m>l` のコストの一部 (例: 係数 α ×
  level m の oct 数) を足し、分割が level 間で入れ子に近づくようにする。各 level の
  均等化は少し崩れるが、粗い level の ghost (細かい level の担当領域の footprint) が
  減る。DICE の円盤 (level 6-7 で ghost が 2-3 倍) や zoom の遷移 level に効く。
  実装場所: `load_balance.f90` の `cmp_new_bound_key_per_level` (コスト `wflag` の計算)。
- 適用 level の限定オプション (`balance_level_min` のようなもの)。効果は限定的。
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
