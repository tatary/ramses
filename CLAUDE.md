# CLAUDE.md — このリポジトリで作業するときの約束

RAMSES (AMR 宇宙流体コード, Fortran 90 + MPI) の岡本研フォーク。
銀河形成シミュレーション (DICE 初期条件, RT, 非平衡化学, 星形成/フィードバック) に使う。

## 最重要ルール: 変更は `patch/takashi/` 以下だけ

- `amr/ hydro/ pm/ poisson/ rt/ io/` は上流 RAMSES のコピー。**直接編集しない**。
- 変更したいファイルは `cp amr/xxx.f90 patch/takashi/` してからコピーを編集する。
  `bin/Makefile` の `VPATH` が `DICE_MHD:DICE:PATCH/$(SOLVER):PATCH:...:../amr` の順に
  探すので、patch 側が自動的に優先される。
- `patch/takashi/to_dice/` (DICE 初期条件用) は `patch/takashi/` より**さらに優先**される。
  `init_part, particle_tree, rho_fine, read_params, init_refine, init_flow_fine,
  init_time, add_list, gadgetreadfile` はここにもあるので、両方を同じ内容に保つ。
- `patch/takashi/mhd/` は `SOLVER=mhd` のときだけ有効。
- `patch/mizuki/` は小野さんの独立した作業コピー。**消さない、こちらから編集しない**。
  差分の取り込みは mizuki → takashi の向きで行う (`diff -w -i -B` で比較)。
- **`patch/takashi/` に `.f90` を置くだけで (git 未追跡でも) ビルドに入る**。使うものは
  `git add`、使わないものは `.f90.mizuki` のように拡張子を変えて外す
  (例: `rt_cooling_module.f90.back`)。
- patch 側のファイルを外した/戻したときは `bin/*.o` の該当オブジェクトを消して
  再コンパイルさせる (make は VPATH 先の古い .o を再利用してしまう)。

## ビルド

```
cd bin && make -j1            # 並列 make はモジュール依存が書かれておらず失敗する
```

- モジュール (`amr_commons`, `rt_parameters` 等) の中身を変えたときは `bin/` の `.o` と `.mod` を
  全部消してからビルドする (依存が書かれていないので古い `.mod` との不整合で落ちる)。

- 本番設定は `bin/Makefile` (RT=1, NIONS=4, NGROUPS=7, NPSCAL=1, DICE=../patch/takashi/to_dice)。
  クラスタ用は `Makefile.xd*` (`make_gnu.sh`, `make_intel.sh`)。
- Mac (Homebrew gfortran) ではリンク前に
  `export LIBRARY_PATH=$(xcrun --show-sdk-path)/usr/lib` が必要。
- 別ディレクトリでビルドするときは Makefile をコピーして `../` を絶対パスに置換する。
- `DEBUG=1` は `-ffpe-trap=...underflow...` を含み RT 初期化が SIGILL で落ちる。
  デバッグ時は `underflow` を外す。
- `LIGHT_MPI_COMM=1` (通信バッファを軽くするフラグ、既定 0) では `reception(icpu,ilevel)` が別の型になる:
  `%igrid/%u/%f` は `%pcomm%igrid/%pcomm%u/%pcomm%f`、粒子の属性は `pcomm%u(property, particle)` と
  添字の順が heavy と逆、送信側は `emission(ilevel)` / `emission_part(ilevel)` に全相手分がまとまる。
  `reception(...)%igrid` や `%up` を触るコードを足すときは必ず `#ifdef LIGHT_MPI_COMM` の両分岐を書き、
  `pcomm` は `if(.not.associated(...))allocate(...)` してから中身を確保する。両方のビルドで出力の bit 一致を
  確認する (DICE: `tests/dice_rt_amr.nml`、grafic: `~/ramses_ks_20230630/ics/ics_L1Mpc_64`、TODO.md 1f)。
- **DICE を有効にしたビルドでは `&INIT_PARAMS` の region 型初期条件は動かない**
  (`to_dice/init_flow_fine.f90` が DICE 専用)。hydro/RT の単体テストは `DICE =` を空にする。

## テスト (patch/takashi/tests/)

| namelist | 内容 | ビルド |
|---|---|---|
| `sedov3d_amr.nml` | Sedov, levelmin 5 / levelmax 7, nremap 2, 4 ランク | RT=0, DICE 空 |
| `stromgren_amr.nml` | Stromgren RT (7 群) を AMR にしたもの | RT=1, DICE 空 |
| `dice_amr.nml` | DICE 矮小銀河 (粒子 194 万), levelmax 8, 星形成あり, RT なし | 本番設定 |
| `dice_rt_amr.nml` | 本番 `xd_rt.nml` の短縮版 (RT, 星の放射, cloudy 冷却, nremap=1) | 本番設定 |

- 回帰確認は「変更前バイナリとの出力バイナリ比較」で行う:
  `git worktree add <tmp> HEAD` で変更前をビルドし、`hydro_*.out`, `part_*.out`,
  `grav_*.out`, `rt_*.out` を `cmp`。分割が変わる変更は `tests/hydro_stats.py`
  (全セルの総和・min/max) で比較する。
- DICE 系は実行ディレクトリの親に `dice/` (dice_dw1.g2), `UVB/`, `cloudy/`, `SED/` (bc03_Chabrier)
  (hm12_ramses.dat) へのシンボリックリンクを置いて使う。データは Google Drive の
  `マイドライブ/codes/ramses_data/` にある (本番 namelist は /mnt/synology1 を指す)。

## namelist の落とし穴

- NIONS=4 ビルドで `rt=.false.` かつ `neq_chem=.false.` だと `rt_init` が落ちる
  (ixHII 未設定)。RT か neq_chem のどちらかを有効にする。
- `rt_err_grad_n` などは `&RT_PARAMS`、`&REFINE_PARAMS` ではない。
- `sf_virial=.true.` は NPSCAL=3 (ivirial1/2) を要求する。
- DICE で `metal=.true.` のとき `zp` は `star.or.sink` のときしか確保されない
  (`eps_star=0` にすると `init_flow_fine` で未確保配列を読む)。
- `initfile`, `uv_file` 等は CHARACTER(80): 長いパスは切れるので相対パスか短い symlink。
- 並列度を変えて restart したいときは Mera.jl の `repartition` で output を書き直す (TODO.md 4 節)。

## コードの地図 (今回の作業に関係する所)

- 分割: `load_balance.f90` (`cmp_new_cpu_map`, `cmp_new_bound_key_per_level`,
  `cmp_cpumap(x,c,nn,ilevel)`, `cmp_cpumap_leaf`), `virtual_boundaries.f90`
  (`authorize_fine`, `build_comm`, 通信ルーチン), `refine_utils.f90` (`make_grid_fine`)。
- level ごとの分割 (balance_per_level) は 2026-10-29 に断念してコードを削除した。経緯と理由は
  `patch/takashi/docs/level_domain_decomposition.md` (記録のみ)。通信の集約・非同期化
  (`make_virtual_fine_dp_nv`, `make_virtual_reverse_dp_start/end`) は残っている。
- 通信: `make_virtual_fine_dp[_nv]` (owner → ghost), `make_virtual_reverse_dp[_nv]`
  (ghost → owner, 加算), `make_virtual_reverse_dp_start/end` (非同期),
  `make_virtual_fine_dp_start/end` (owner → ghost の非同期, RT subcycle の `rtuold` で使用)。
  reverse の start/end の間で level ilevel の仮想セルを書き換えてはいけない。forward は start 時点の
  owner の値が送られる (送信バッファに詰めるので start 後の書き換えは可) が、end までは仮想セルが古い。
  `cooling_fine_part(ilevel,1/2)` は emission grid (他ランクが ghost として持つ自分の grid) を先に
  冷却するための分割版 (`cooling_fine` と同じ結果)。`star_formation` の `uold(1:4)` 交換は
  `amr_step` が slot 1 で start し、`star_formation` が end して ghost を自分で速度に変換する。
- ks (杉村さん版) からの移植 (`merge_ks` ブランチ、`docs/merge_ks.md`): 星形成は `&SF_PARAMS` の `sf_fg=.true.` で
  `star_formation_fg.f90` (Pop II 局所/SFC、Pop III、非局所枯渇) に切り替わる。このとき粒子 family 3-8 は
  `FAM_POP3/SN/DEAD/BH/SFC/PSC` の意味になる (`pm_commons`; 既定では上流の cloud/debris/other)。補助: `my_utils.f90`
  (ufront)、`sbs_imf.f90`、`output_analysis.f90` (ugvtk: SFC/Pop III/指定星の周りの入れ子 VTK、`&OUTPUT_PARAMS` の `ugvtk=.true.` で有効)、`sf_log.f90` (星形成/SN ログの共通 writer)。
  テストは `tests/dice_rt_fg.nml`。フィードバック側: `feedback.f90` に Pop III SN と SFC 寿命、`&FEEDBACK_PARAMS` の
  `SN_reservoir` (SN 熱を冷却半径の球へ、`sn_*.f90` 3 本 + amr_step/refine_utils/load_balance/init_hydro/output_hydro のフック)、
  `fcool_stop` (冷却停止タイマー、NPSCAL を 1 つ余計に取る)。sink/BH: `sf_fg` のとき `sink_particle_fg.f90` (ks の
  sink、全ルーチン `_fg` 付き、cloud 粒子は family 9 `FAM_CLOUD_FG` で `is_cloud` が sf_fg を見る) に amr_step/init_part/
  particle_tree/output_amr/newdt_fine/force_fine/read_params が振り分ける。`&SINK_PARAMS` の `bhpop3`, `BH_model`,
  `sink_mass_min_pop3BH` 等; BH の放射は `rt_spectra` の `set_BH_SED`/`BH_RT_feedback`。
- やり残し・既知のバグ: `patch/takashi/TODO.md`。

## 作業の進め方

- 物理結果を変えない変更 (通信順序の入れ替え等) は出力の bit 一致で確認する。
- 大きな変更は設計メモを `patch/takashi/docs/` に残す (日本語でよい)。やり残しは `patch/takashi/TODO.md` に集約する。
- commit メッセージ・PR は日本語でも英語でも可。commit は頼まれたときだけ。
