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
  (例: `flag_utils.f90.mizuki`, バグあり。TODO.md 参照)。
- patch 側のファイルを外した/戻したときは `bin/*.o` の該当オブジェクトを消して
  再コンパイルさせる (make は VPATH 先の古い .o を再利用してしまう)。

## ビルド

```
cd bin && make -j1            # 並列 make はモジュール依存が書かれておらず失敗する
```

- 本番設定は `bin/Makefile` (RT=1, NIONS=4, NGROUPS=7, NPSCAL=1, DICE=../patch/takashi/to_dice)。
  クラスタ用は `Makefile.xd*` (`make_gnu.sh`, `make_intel.sh`)。
- Mac (Homebrew gfortran) ではリンク前に
  `export LIBRARY_PATH=$(xcrun --show-sdk-path)/usr/lib` が必要。
- 別ディレクトリでビルドするときは Makefile をコピーして `../` を絶対パスに置換する。
- `DEBUG=1` は `-ffpe-trap=...underflow...` を含み RT 初期化が SIGILL で落ちる。
  デバッグ時は `underflow` を外す。
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
- `balance_per_level=.true.` を `&RUN_PARAMS` に足すと level ごとの分割になる。
  ghost を抑えるオプションは `balance_level_min`, `balance_tolerance` (docs 4c 節)。

## namelist の落とし穴

- NIONS=4 ビルドで `rt=.false.` かつ `neq_chem=.false.` だと `rt_init` が落ちる
  (ixHII 未設定)。RT か neq_chem のどちらかを有効にする。
- `rt_err_grad_n` などは `&RT_PARAMS`、`&REFINE_PARAMS` ではない。
- `sf_virial=.true.` は NPSCAL=3 (ivirial1/2) を要求する。
- DICE で `metal=.true.` のとき `zp` は `star.or.sink` のときしか確保されない
  (`eps_star=0` にすると `init_flow_fine` で未確保配列を読む)。
- `initfile`, `uv_file` 等は CHARACTER(80): 長いパスは切れるので相対パスか短い symlink。

## コードの地図 (今回の作業に関係する所)

- 分割: `load_balance.f90` (`cmp_new_cpu_map`, `cmp_new_bound_key_per_level`,
  `cmp_cpumap(x,c,nn,ilevel)`, `cmp_cpumap_leaf`), `virtual_boundaries.f90`
  (`authorize_fine`, `build_comm`, 通信ルーチン), `refine_utils.f90` (`make_grid_fine`)。
- level ごとの分割の設計・検証結果: `patch/takashi/docs/level_domain_decomposition.md`。
- 通信: `make_virtual_fine_dp[_nv]` (owner → ghost), `make_virtual_reverse_dp[_nv]`
  (ghost → owner, 加算), `make_virtual_reverse_dp_start/end` (非同期)。
  reverse の start/end の間で level ilevel の仮想セルを書き換えてはいけない。
- やり残し・既知のバグ: `patch/takashi/TODO.md`。

## 作業の進め方

- 物理結果を変えない変更 (通信順序の入れ替え等) は出力の bit 一致で確認する。
- 大きな変更は設計メモを `patch/takashi/docs/` に残す (日本語でよい)。やり残しは `patch/takashi/TODO.md` に集約する。
- commit メッセージ・PR は日本語でも英語でも可。commit は頼まれたときだけ。
