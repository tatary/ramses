# patch/ks (Ramses-FG) の物理・最適化を patch/takashi に取り込む計画

作成 2026-11-12。ブランチ `merge_ks` (`new_dec` から分岐)。取り込み元は
`/home/takashi/work7/ramses_test/Ramses-FG` の `dev_to` ブランチ、コミット 40a68ca (固定)。

## 前提 (調べて分かったこと)

- 両リポジトリの上流 (`amr/ hydro/ pm/ poisson/ rt/`) は 1 ファイルも差がない。両 patch は同じ土台への
  独立した差分なので、上流を共通祖先にした 3-way マージ (`git merge-file takashi 上流 ks`) が使える。
  git 履歴は共有していない (FG は 2019 年のスナップショットから 213 コミット)。
- ks の `.f90` は 86 本: 上流にもある 62 本 (うち takashi も触っている 41 本) と新規 24 本
  (`mechanical_fine*`, `sn_*`, `fof*`, `sbs_imf`, `output_analysis`, `output_vtk`, `rt_background`,
  `spread_fine`, `nonlocal_feedback`, `feedback_flag`, `add_rho`, `posix_fs`, `my_utils` など)。
- ks の差分は物理に集中: `star_formation` (上流比 7159 行)、`sink_particle` (3299)、`feedback` (2952)、
  `rt_cooling_module` (2234)、`rt_spectra` (2056)、`amr_step` (1568)、`coolrates_module` (1261)、
  `load_balance` (966)。takashi の差分は通信 (`virtual_boundaries` 1311、`update_time` 560、
  `multigrid_fine_commons` 448)、DICE (`init_*`, `rho_fine`, `particle_tree`)、OpenMP、`cg_auto`。
  ks には通信の集約・非同期化、OpenMP、`LIGHT_MPI_COMM` は入っていない。
- ks のビルド: `makefiles/Makefile_0X.xd`、NVAR=12、NGROUPS=5 (4 + LW)、NIONS=4、QUADHILBERT、DICE なし。
  ks の namelist 例は `namelist/test_problem/smbh/ngridmax_fire3.nml`、本番ジョブは
  `jobfiles/xd/smbh/ngridmax_intel.sh` (sed で上書きする項目あり)。
- **整形**: 両 patch とも `findent -i2 -Rr` で整形済み、上流は未整形。そのまま 3-way すると衝突 ~900 か所、
  3 者とも findent を通すと 165 か所 (9,400 行) に減る (`init_part` 80→10、`rt_spectra` 77→3)。
  多い順: `rt_cooling_module` 29、`feedback` 16、`amr_step` 14、`star_formation` 14、`coolrates_module` 13、
  `virtual_boundaries` 11、`init_part` 10。試算スクリプトは scratchpad `merge_dry2/` (findent は同 venv)。

## 方針

1. **metal cooling table (cloudy) と UV background は takashi 版を優先。** ks の `rt_cooling_module` には
   cloudy テーブル (`init_cloudy_metal_cooling`, `rt_cmp_cloudy_metals`, `ss_factor_rahmati`) が無く、
   別実装 `rt_cmp_metals_cloudy` を使っている。cooling まわりで ks から持ってくるのは H2/LW 化学、dust、
   2 次電離、Pop III 放射源、subcycle の最適化 (xHI の 10% 則を MIN 形に; NOTES_rt_cooling_highT_subcycle.md §6)。
2. **分子 (H2) 冷却は ks 版の高密度対応を取り込み、namelist で有効化。既定は上流のまま (bit 一致)。**
   → 2026-11-12 に移植済み (下の「cooling の移植」)。
3. LW 群の添字 `LWgroup` は namelist で明示指定 (既定 -1 = LW なし)。ks の「`isH2Katz` なら `LWgroup=ixHI`」
   という自動設定は持ち込まない。`xd_rt.nml` の 7 群構成では群 4 (11.2–13.6 eV) が LW なので `LWgroup=4`。
   → 2026-11-12 に移植済み (下の「LW 群自己遮蔽の移植」)。
4. ks の機能は全部 namelist スイッチで既定オフ。各段階で `tests/` の 4 本が変更前バイナリと bit 一致。
5. ks 用の回帰テストを用意する (ks の設定を takashi の Makefile で組む: NVAR→NPSCAL、NGROUPS=5、QUADHILBERT、
   DICE 空; 小さな FG 初期条件で FG バイナリと bit 比較)。
   2026-11-12 に ks の本番 `ngridmax_intel.sh` に対応する 3 ファイルを作った (下の「ks 本番相当の設定」)。
   FG 初期条件での比較ランは未実施。
6. `~`, `.orig`, `.back.f90`, `.nozaki.f90`, `backup_*`, `stdout*` は持ち込まない。
7. 物理定数が食い違ったら精度の高い方 (`amr/constants.f90`, CODATA 2018) を使う。

## 手順 (段階ごとにビルドと bit 一致を確認)

1. モジュールと namelist: `amr_parameters`, `amr_commons`, `pm_commons`, `rt_parameters`, `hydro_parameters`,
   `read_params`, `read_hydro_params`, `rt_init`、新規 `mechanical_commons`, `fof_commons`。
2. ks だけが触っているファイル (上流ファイル約 20 本 + 新規 24 本) は衝突なしでコピー。
3. 基盤: `amr_step`, `virtual_boundaries`, `update_time`, `init_part`, `load_balance`, `output_amr`。
   `to_dice/` の 8 ファイルは ks が全部触っているので両方に同じ変更を入れる。
4. 物理: `star_formation`, `feedback` + `mechanical_fine` + `sn_*`, `sink_particle`, Pop III。
5. 冷却と RT: `rt_cooling_module` (最難関: 同じ cloudy コードから両側が別々に進化、ks の高温 subcycle と
   takashi の OpenMP が同じ箇所)、`coolrates_module`, `cooling_fine`, `rt_flux_module`, `rt_spectra`。

## cooling の移植 (2026-11-12、方針 2)

`coolrates_module.f90` / `amr_parameters.f90` / `read_hydro_params.f90` / `rt_cooling_module.f90` (呼び出し 1 行)。
`&COOLING_PARAMS` に次を追加 (ks と同じ名前、既定は上流の挙動):

| 変数 | 既定 | 内容 |
|---|---|---|
| `H2coolfunc` | `''` | `''`/`'Hollenbach'`: HM79+HC12 の低密度極限 (上流のテーブル、LTE 頭打ちなし)。`'Galli'`: Galli & Palla 1998 (H 衝突、HM79 の LTE 極限)。`'Glover'`: Glover & Abel 2008 (H, H2, He, H+, e-) + Glover 2015 LTE。ks 本番は `'Glover'` |
| `H2_fesc` | `.false.` | H2 線冷却に脱出確率 (Fukushima, Omukai & Hosokawa 2018; Jeans 長の半分 × n_H2 の柱密度) を掛ける。ks は常時オン |
| `H2_Hm_cool` | `.false.` | H⁻ 放射付着冷却 (H + e → H⁻ + hν)。ks は常時オン (`no_H2Hm_cool=.false.`) |
| `H2_chem_cool` | `.false.` | H2 生成・解離の化学的加熱・冷却 (H⁻ 経路、三体、衝突解離 `comp_Beta_H2coll`)。`H2_Hm_cool` を自動で立てる。ks 本番はオフ |

- `compCoolrate(T, ne, nH, mu, nN, nI, dcooldT)` に `nH, mu` を追加 (H2 の新オプションだけが使う)。
- 既定経路は式の順序を変えていない (`cr_H2HI`, `cr_H2H2` の加算順そのまま、脱出確率の乗算と H⁻ 項の加算は
  スイッチが立ったときだけ) ので bit 一致。
- 新オプションは `T < T_H2chem_max = 1e6 K` でだけ評価し、H⁻ のフィットは `min(T, T_H2chem_fitmax = 1e5 K)` で
  評価する (ks の TO MODIFIED の扱いを踏襲)。
- 定数は `constants` モジュール (kB, mH, eV2erg, factG_in_cgs, twopi) に統一。ks の手書き値
  (1.380662e-16, 1.6022e-12, 6.674e-8) は使わない。
- `'Galli'` の `LambLTE/nHI` を使う形は nHI=0 で割り算になるので、代数的に同じ `nH2 L_LTE nHI Λ0/(nHI Λ0 + L_LTE)` に。
- 持ち込まなかったもの: `fclump_H2dust` (H2 の dust 生成率の clumping、化学側)、`T2_lfloor`、`cool_step_dbg` の
  デバッグ出力、`no_H2Hm_cool` (既定の向きが逆なので `H2_Hm_cool` に置き換え)、LW 群の自己遮蔽
  (`isH2Katz`, `LW_ss_*`; RT 群構成に関わるので別段階)。

## ks 本番相当の設定 (2026-11-12、方針 5)

ks の `jobfiles/xd/smbh/ngridmax_intel.sh` (ベース `namelist/test_problem/smbh/ngridmax_fire3.nml` に sed で
上書き) に対応するもの:

- `bin/Makefile.xd_ks_intel`: `Makefile.xd_intel` との差は DICE 空、NGROUPS=5、NPSCAL=3 (NVAR=12: 6 metal、
  7-10 ions、11 zoom scalar、12 は ks の cooling flag 用の予約)、`QUADHILBERT=1` (フラグを追加)、EXEC=ramses_ks。
  bin/ と .o が混ざるので `bin_ks/` にコピーして組む (相対パスはそのまま通る)。gfortran でビルド確認済み。
- `namelist/ks/ngridmax_fire3.nml`: ks の値 (sed 上書き後) を takashi の namelist 群に翻訳。takashi に無い
  項目は `!KS` でコメント (106 行; 移植したら `!KS` を外す)。意味が違う項目は `!TO` で注記: `m_star`
  (ks は Pop II 粒子の最小質量 [Msun])、`n_star` (ks は負で Jeans 長閾値)、`pressure_fix_cond` (ks は 3)、
  `walltime_hrs` (ks は `maxruntime`)、`sink`/`create_sinks` (ks は bhpop3; 上流の clump ベースと別物なのでオフ)、
  `PEH` (takashi は `iPEH_group`)、`rt_c_fraction` (ks は level ごと)。RT 群は ks の 5 群 (LW, HI, HeI, HeII, FUV)
  で、断面積は `xd_rt.nml` の同じ帯の値。cloudy テーブルは方針 1 で takashi のもの。
  1 ランクで読み込み確認済み (初期条件ファイルが無いところまで進む)。
- `jobfiles/xd/ks/ngridmax_intel.sh`: ks と同じ SBATCH (4 ノード 448 ランク、P-large-a、24 h)、restart 自動判定、
  `restart_files` の退避。binary は `bin_ks/ramses_ks3d`、初期条件・SED・cloudy は symlink (`ICS`, `SED`, `cloudy/`) で
  相対パスにする (CHARACTER(80) 対策)。flat ラン; hybrid にするなら `run_hybrid_intel.rb` の設定を足す。
- ks の Makefile は 4 つとも `-DQUADHILBERT` を DEFINES に直書き (スイッチなし)。

## ks の 2 層 domain decomposition (two_tier_balance) の移植性 (2026-11-12 調査)

状態: ks では実装済み (2026-11-05, `TO ADDED/MODIFIED 2026-11-05/06`)、クラスタで検証 PASS (11/10)、本番投入済み
(dfaa7e7, `two_tier_balance=.true., two_tier_level=13`)。設計は Dropbox 側 (patch/takashi) で断念した per-level 分割の
後継で、per-level のバグ 4b-2 (孤児 grid) は `test_flag` の修正で潰している。

触っているファイルとサブルーチン (findent で揃えた上で上流との差をサブルーチン単位で比較):

| ファイル | two-tier の変更 | takashi 側との重なり |
|---|---|---|
| `amr_parameters`, `read_params` (+`to_dice/read_params`) | namelist 4 変数 (`two_tier_balance/level/fine_frac/test_shift`) | 追記のみ |
| `amr_commons` | `bound_key_f/2_f`, `two_tier_ls/2`, footprint 配列、module 関数 4 つ (~60 行) | 追記のみ |
| `load_balance` | `cmp_new_cpu_map` 冒頭の分岐 → 新規 `cmp_new_cpu_map_two_tier`, `two_tier_cut`, `two_tier_*_foot`, `cmp_cpumap_lev`, `cmp_cpumap_pos`, `two_tier_report`, `two_tier_cpu_hull` (新規 ~560 行) と `load_balance` の末尾 (bound_key_f の更新) | takashi の `cmp_new_cpu_map` (cost_weight_power, cost_time_weighting) と同じルーチンだが別の場所。**この段階では持ち込まないもの**: `cost_weighting=.false.` 時の `niter_cost=2*(ilevel-levelmin)` (上流と違う既定、以後も不要)、angular ordering (不要)、SN reservoir の転送 (`sn_pending` のフック; `sn_energy_reservoir` が無いと組めないので手順 4 のフィードバック移植で一緒に入れる。`load_balance`, `defrag`, `refine_fine`, `kill_grid`, `make_grid_fine` の 5 か所) |
| `virtual_boundaries` | `authorize_fine` (two-tier 分岐 ~100 行、上流部分は `if(.not.two_tier_balance)` に残る)、`build_comm` の index-0 検査、新規 `two_tier_authorized`, `check_tree_consistency`, `debug_missing_owner`, `get_mem` | takashi は `authorize_fine` を触っていない。`build_comm` は両方が触るが ks 側は検査の追加だけ (ks の make_virtual_* の差分は LIGHT_MPI_COMM の削除で実質なし) |
| `refine_utils` | `make_grid_coarse` / `make_grid_fine` の `cmp_cpumap` → `cmp_cpumap_lev(…,2)` / `(…,ilevel+1)` の 2 行 | takashi は上流のまま → 2 行だけの patch を新設 |
| `flag_utils` | `test_flag` の修正 (子の「refined」を所有者から flag2 で配る、14 行、ks では常時有効) | takashi の `userflag_fine` とは別ルーチン |
| `init_amr` | 初期値 (`two_tier_ls=nlevelmax`)、level-1 grid の所有者 `cmp_cpumap_lev(…,1)`、restart の末尾レコード `read_two_tier_record`、`two_tier_bcast_foot` | takashi の init_amr (DICE) と同じルーチンだが別の場所 |
| `output_amr` | ordering レコード `'hilbert2t'`、末尾レコード (`two_tier_ls`, footprint)、`info_twotier_XXXXX.txt`、DOMAIN 表の hull | takashi の output_amr (dump_restart 等) とは別の場所 |
| 位置→ランク | `sink_particle` (2), `rt_spectra` (BH 光子), `amr_step` (dbg), `my_utils` → `cmp_cpumap_pos/_lev` | ks 物理の側。takashi では `init_part`/`to_dice/init_part` (粒子初期配置) と `init_amr` が `cmp_cpumap` を呼ぶが、粒子配置は粗い層のままでよい (ks も同じ) |

注意点:
- 出力形式が変わる (`hilbert2t` + 末尾レコード)。Mera の `repartition`/読み込みは非対応なので、two-tier の output からは
  ランク数変更ができない (ks は単一分割の restart を別に保持)。既定 `.false.` では出力は不変。
- `test_flag` の修正は ks では常時有効 (16 ランクの試験では bit 一致、本番では境界例で buffer oct が増える方向)。
  takashi でも常時有効にした (2026-11-13): 単一分割でも踏める孤児 grid のバグなので two-tier 専用にしない。
  既定テスト (dice_rt_amr heavy/LIGHT、dice_amr) は常時有効でも変更前バイナリと bit 一致。
- `-DTWO_TIER_CHECK` の集団診断 (`check_tree_consistency`, `debug_missing_owner`, `check_particle_tree`) も一緒に持ってくる。
- takashi の DICE 本番 (levelmin 6 / levelmax 15) は ks の zoom (levelmin 7 / levelmax 21、grid の 97% が粗い level) と
  構造が違い、balance 後の grids-work max/avg は 1.10 (TODO 2d-3)。効くかどうかは別問題。

### Mera.jl の two-tier / QUADHILBERT 対応 (後回し、2026-11-12 の判断)

- two-tier 出力の読み書き: 読み側は実質そのまま (ordering は "bisection" かどうかしか見ない、所有者はファイルから、末尾レコードは警告のみ)。
  書き側は `'hilbert2t'` を `'hilbert'` に戻して単一分割で書くだけ (RAMSES が restart 後の最初の balance で 2 層に組み直す)。
- QUADHILBERT (ks 本番は全部): `bound_key` が 16 バイトで今は error。levelmax 21 では key ~2^63 で Float64 の仮数 53 bit を超えるので、
  key と区切りを `UInt128` で計算し、ファイルとは IEEE binary128 ↔ UInt128 の変換 (`hilbert.jl`, `decomposition.jl`, `amr_file.jl`)。
- yt (4.3.1) は two-tier 出力をそのまま読める (ordering レコードは読み捨て、末尾レコードは読まない、info の DOMAIN 表は hull で上位集合)。

## LW 群自己遮蔽の移植 (2026-11-12、方針 3)

ks の Katz 型 H2/LW モデル (`isH2Katz`) は、上流 (Nickerson+2018) の H2 モデルと同じ系統 (LW 光解離断面積、
6.94 の pumping 加熱、Burton 1990 の de-excitation、dust 吸収) に、遮蔽係数を掛けたもの。こちらは上流の
モデルを土台に、遮蔽係数だけを足した。`&RT_PARAMS` (既定は全部オフ = bit 一致):

| 変数 | 既定 | ks の対応 | 内容 |
|---|---|---|---|
| `LWgroup` | -1 | `isH2Katz` で `ixHI`(=1) に固定 | LW 帯の光子群。`xd_rt.nml` では 4 |
| `LW_ss_model` | `'none'` | `LW_ss_wgh=.true.` → `'WGH19'` | H2 自己遮蔽: `'WGH19'` Wolcott-Green & Haiman 2019、`'DB96'` Draine & Bertoldi 1996 (w=0.035)、`'G09'` Gnedin+2009 形 (w=0.2) |
| `LW_ss_dust` | `.true.` | 常時 | dust 遮蔽 exp(-4e-21 Z N_H) (Gnedin+2009) も掛ける |
| `LW_abs_thin` | `.true.` | 同名 | `.true.`: 遮蔽は光解離率だけ。`.false.`: LW 光子の H2 吸収にも |
| `LW_abs_atom` | `.false.` | 同名 (本番 .true.) | LW 光子の HI 吸収 (5.23e-25 cm²) |
| `LW_ss_from_incident_surface_yH2` | `.false.` | 同名 (本番 .true.) | 遮蔽なしの平衡 xH2 = cr/de で遮蔽係数を再計算 (光解離率にのみ使用) |

- 実装: `rt_solve_cooling` に optional 引数 `dx_SS` (遮蔽長 [cm]) を追加し、`cooling_fine` がセル幅 `scale_l*dx_loc` を渡す
  (`rt_evol_single_cell` は渡さない → 遮蔽なし)。`cool_step` で `f_shd = LW_shield_factor(nHI, nH2, nH, TK, Z, dx)` を
  現在の xH2 から計算し、(1) H2 光解離率 `SUM(signc(:,ixHI)*dNp*fshg)` (`fshg(LWgroup)=f_shd`)、(2) LW pumping 加熱の
  同じ和、(3) `LW_abs_thin=.false.` なら `phAbs(LWgroup)` の H2 吸収、に掛ける。`LW_abs_atom` は `phAbs(LWgroup)` に
  `nHI*5.23e-25*c` を足す。incident-surface オプションは H2 更新の直前で `f_shd` を再計算して (1) に使う
  (ks では save 変数経由で次の呼び出しの吸収・加熱にも効いていたが、ここでは同じ呼び出し内の光解離率だけ)。
- 遮蔽関数 `comp_Sd`, `comp_SH2_G09`, `comp_SH2_DB96`, `comp_SH2_WGH` は `coolrates_module` に (ks と同式、定数は constants)。
- 持ち込まなかったもの: ks の薄い極限の断面積 `Dlw/Flw` (こちらは `group_csn(LWgroup, ixHI)`)、`ssh2`=400 の扱いは上流のまま、
  `SS_LVG` (LVG 近似の遮蔽長; ks 本番は .false.)、`H2chem_hot` (高温で H2 化学を飛ばす; ks でも作業ツリーから外してある)。
- 検証: `dice_rt_amr.nml` (既定) と `dice_amr.nml` が変更前と bit 一致 (40 + 32 ファイル)。`tests/dice_rt_amr.nml` は
  `star_as_nbody=.true.` で光子が 0 なので LW の効果は出ない。`star_as_nbody=.false.` にして ks 本番相当
  (`LWgroup=4, LW_ss_model='WGH19', LW_abs_atom=.true., LW_ss_from_incident_surface_yH2=.true.`) で走らせると完走し、
  hydro/rt の 16 ファイル中 8 が変わる (圧力の総和 +7e-6、LW 群の光子数 -4e-5 の相対変化; 6 step・levelmax 8 では H2 が少ない)。

## ks の H2 化学モデル `isH2Katz` の移植 (2026-11-12)

ks の計算をこちらで回すには、cooling だけでなく H2 の生成・解離率そのものが ks と同じでなければならない
(始原ガスでは H⁻ 経路が電離度に比例する点が効く)。上流 (Nickerson+2018) のテーブル方式を既定のまま残し、
`&RT_PARAMS` の `isH2Katz=.true.` (ks と同名) で次を ks の形に切り替える:

| 項目 | 上流 (既定) | `isH2Katz` |
|---|---|---|
| dust 表面生成 | 9e-17 √T₂/f(T₂) Z f_dust (HM79) | 3.5e-17 Z `fclump_H2dust` (1-xHII) (Katz+2017) |
| 気相生成 (H⁻) | 8e-19 (T/1000)^0.88 n_e | Glover+2010 の k1, k2, k5, k13 で平衡 H⁻: k1 k2 x_e x_HI/(k2 x_HI + k5 x_HII + k13 x_e)。フィットは min(T, 1e5 K) で評価 |
| 三体 | テーブル (Forrey 2013, Palla+1983) | 同じ式を直接評価 |
| 衝突解離 | テーブル (Dove & Mandy 1986, Martin+1998) | `comp_Beta_H2coll`: Glover & Abel 2008 の e⁻, H, H2, He, H⁺ 衝突、臨界密度で LTE/非 LTE を混合 |
| T ≥ 1e6 K | そのまま | 生成 0、解離は高温極限の式、光解離を飛ばす (`H2chem_hot`) |
| H2 電離閾値 | 13.6 eV (`ionEv_HI`) | 15.2 eV |
| H2 断面積 (SED から群の csn を作るとき) | Sternberg 2014 の 2.1e-19 (LW) + Abel 1997 | Baczynski+2015 Table 1 (LW は 2.47e-18 = Katz の Dlw/Flw) |

- `comp_Alpha_H2_Katz` を `coolrates_module` に追加 (ks の `comp_Alpha_H2` の isH2Katz 分岐; 定数は constants)。
  `cool_step` では `alpha(ixHI) = comp_Alpha_H2_Katz(...) * nH` として上流の `cr = alpha*xHI` の規約に合わせる。
- `fclump_H2dust` は `&COOLING_PARAMS` (ks と同じ群)。
- ks で `isH2Katz` が同時に切り替えていた LW の扱い (`LWgroup=ixHI` 固定、`sigUVthin=Dlw/Flw`) は、こちらでは
  `LWgroup` と群の csn (SED から Baczynski の 2.47e-18 で更新される) で表現する。`sedprops_update=1` が前提。
- 検証: 既定で `dice_rt_amr` が bit 一致 (40 ファイル)。`isH2Katz` + Glover + fesc + H⁻ + LW 遮蔽 + 放射する星で完走、
  圧力総和 -0.13%、xHI +0.18%、Coolstats の平均ループ 1.035 / 最大 184 (比較対象 1.036 / 225)。

## cooling の subcycle (2026-11-12 調査)

ks の修正 A〜H を 1 つずつ突き合わせた結果。ks の計算をこちらで回すには、ks の領域 (SN 加熱セル、始原ガス、電離フロント) で
効いた修正がこちらにも要る。

| ks | 内容 | こちら |
|---|---|---|
| A | H⁻/H2 フィットを min(T, 1e5 K) で評価 | `T_H2chem_fitmax` として H2 冷却 (`chemical_and_Hm_cool`, `comp_Beta_H2coll`) と `comp_Alpha_H2_Katz` に入れた |
| B | T > 1e6 K で H2/H⁻ 冷却を飛ばす | `T_H2chem_max` で同じ |
| C | `dxL_H2`/`dxL_HM` の微分のタイポ | 新規実装なので該当なし |
| D | ddt アンダーフローのフォールバック | ks でも未着手 |
| E | T > 1e6 K で H2 化学も飛ばす (`H2chem_hot`) | `isH2Katz` に入れた (ks の作業ツリーでは外してあるが、ノートの結論は「物理的に妥当」) |
| F | xHI の 10% 則を MIN(\|Δx_HI\|/x_HI, \|Δn_e\|/n_e) に | `relax_x_by_ne` (commit 2a712c51) が H と He の全種に同じ形で既にある。既定 `.false.`。ks 用 namelist では `.true.` |
| G | n_e 判定の単位 | こちらは `neInit + x_FM*nH` で元から正しい |
| H | xHII 更新の生成項 `cr*(1-xHII)` → `cr*xHI` | こちらの更新は元から `cr=(beta*ne+photoRate)*dxHI` |
| — | Glover_rate の e⁻ ブロック 2 件 | 移植済み |
| — | Coolstats の MPI 集約 (INTEGER8, SUM) と T2 判定の NaN 対策 | 既に入っている (TO ADDED) |
| — | ループ上限 | こちら 1e5 (上流)、ks は 1e6。ks の修正後は最大 256 なので 1e5 のまま |

こちらの DICE 本番 (flat-334557) の Coolstats は平均ループ 1.0035、最大 660、cooling は step 104 s 中 9 s。
`relax_x_by_ne=.true.` の効果はこちらでは未検証 (ks の検証では物理不変)。

## 光電加熱の統一 (2026-11-12): こちらの `iPEH_group` 方式に

ks の `PEH` (Kim & Ostriker 2017 形: 2.8e-26 G0 n_H Z、効率は定数 0.0215 相当、G0 は 5.4-13.6 eV の全群、T ≤ T_sputter で
dust あり、光子は消費しない) は持ち込まず、こちらの方式 (Bakes & Tielens 1994 の ε(G0,T,n_e)、1.3e-24 ε G0 n_H Z (1-xHII)、
吸収分を群から除く) に統一する。LW 帯 (11.2-13.6 eV) も G0 に数えるのが Habing/Bakes & Tielens の定義どおりなので、
`PEH_all_FUV=.true.` (既定 .false.) で `5.4 < group_egy <= 13.6 eV` の全群の和を G0 にし、吸収も同じ断面積で各群に掛ける
拡張を入れた。ks 用 namelist は `iPEH_group=5, PEH_all_FUV=.true.`。ks との残差は ε が定数でないこと (弱い放射場で
2.3 倍、強い帯電で小さい) だけ。

## 2 次電離 `sec_ion` の移植 (2026-11-12)

ks には 2 つの実装があり、Park さんのテーブル版 (`fifh`、UV 背景用、`uv_file` の特別な形式) は cool_step では
コメントアウトされていて使われていない。生きているのは KS 版で、Ricotti, Gnedin & Shull 2002 のフィット
(`secion_HI2nd`, `secion_HeI2nd`, `secion_heat`: 1 次電子のエネルギー E0 と電離度 x_e = n_e/(n_H+n_He) の関数) を
HI と HeI の光電離に適用する。こちらにはこれを `&RT_PARAMS` の `sec_ion` (既定 .false.) として入れた。

- `Eprim_avg(群, イオン) = sigec*egy/signc - ionEv` (1 次電子の平均エネルギー [eV]) を `updateRTGroups_CoolConstants` で計算。
- 加熱: HI と HeI の光電離の分は `PHrate` の代わりに `signc * Eh(E0, x_e)` (熱になる割合)。H2 と HeII は `PHrate` のまま。
- HI の電離率に `signc(HI)*Phi_HI(E0_HI) + signc(HeI)*Phi_HI(E0_HeI)*n_HeI/n_HI`、HeI も同様に 2 次電離を加算。
- ks は `igroup >= Xgroup` (X 線群) だけに適用していた (加熱側は `.gt.`、電離側は `.ge.` と不揃い)。こちらは全群に適用する:
  フィット自体が E0 < 28 eV で 2 次電離 0、E0 < 11 eV で全部熱、を返すので、低エネルギー群では従来式と一致し
  (signc*E0 = sigec*egy - signc*ionEv = PHrate)、群の添字に依存しない。
- ks の `stop` (x_e の範囲外) は x_e のクランプに置き換え、関数は pure に。
- 検証: 既定で bit 一致 (40 ファイル)。`sec_ion=.true.` + 放射する星で完走、圧力総和 -2e-5 (星の SED は 65.7 eV の群まで
  なので 1 次電子が 28 eV を超えるのは HeII 群の一部だけ)。

## 2 層分割 (two_tier_balance) の移植 (2026-11-13)

上の「移植性」の表のとおりに入れた。namelist (`&RUN_PARAMS`): `two_tier_balance` (既定 .false.)、`two_tier_level`
(Ls、0 で自動)、`two_tier_fine_frac`、`two_tier_test_shift`。`ordering='bisection'` とは併用不可 (read_params で停止)。

- `amr_parameters`, `read_params` (+`to_dice`), `amr_commons` (`bound_key_f/2_f`, `two_tier_ls/2`, footprint 配列、
  `contains` に補助関数 5 つ): ks のまま。
- `load_balance.f90`: `load_balance` の末尾 (bound_key_f と footprint の更新、`two_tier_report`、STATUS 行)、
  `cmp_new_cpu_map` 冒頭の分岐、ks の新規ルーチン 14 本 (`cmp_new_cpu_map_two_tier`, `two_tier_cut`, footprint 処理、
  `cmp_cpumap_lev`, `cmp_cpumap_pos`, `two_tier_report`, `check_particle_tree` など) を末尾に追加。`check_cloud_count` (sink) は
  入れていない。ks の `include 'mpif.h'` は `use mpi_mod` に (gfortran は同一ファイル内の型違いの暗黙呼び出しを拒否する)。
  key (`real(qdp)`) の通信は `MPI_BYTE` なので QUADHILBERT でもそのまま (ks 設定のビルドで確認)。
- `virtual_boundaries.f90`: `authorize_fine` に 2 層の分岐 (現在と新しい分割の両方)、`build_comm` の末尾に
  「他ランクが ghost として持つ自分の grid が手元に無い」検査 (heavy/LIGHT 両方)、`-DTWO_TIER_CHECK` で集団診断
  (`check_tree_consistency`, `two_tier_authorized`, `debug_missing_owner`; heavy 構造のみなので `#ifndef LIGHT_MPI_COMM`)。
- `refine_utils.f90` を patch に新設 (上流 + `cmp_cpumap_lev` の 2 行)。
- `flag_utils.f90` `test_flag`: 「子が他ランクで refine されている」を所有者から flag2 で配る修正。ks と同じく常時有効
  (当初 `two_tier_balance` でゲートしていたが、単一分割でも踏めるバグなので外した; 既定テストは bit 一致のまま)。
  コストは level ごとに `make_virtual_fine_int(flag2, ilevel+1)` が 1 回増えるだけ。
- `init_amr.f90`: 初期値 (`bound_key_f=bound_key`, `Ls=nlevelmax`)、level-1 grid の所有者 `cmp_cpumap_lev(…,1)`、restart の
  `'2t'` 印の検査と末尾レコード (`read_two_tier_record`; 2 系統の reader の両方)、`two_tier_bcast_foot`。
- `output_amr.f90`: amr ファイルの ordering レコードを `'hilbert2t'` に、末尾レコード (Ls, bound_key_f, footprint)、
  `info_twotier_XXXXX.txt`、info の DOMAIN 表は hull。
- ks から持ち込まなかったもの: `cost_weighting=.false.` 時の線形重み、SN reservoir の転送 (手順 4 で)、
  `sink_particle`/`rt_spectra`/`my_utils` の `cmp_cpumap_pos` 呼び出し (ks 物理側)、`get_mem`。

検証 (dice_rt_amr, 4 ランク, levelmin 6 / levelmax 8, 6 step, nremap=1):
- 既定 (`.false.`): heavy と LIGHT のビルドが変更前と bit 一致 (40 ファイルずつ)、`dice_amr` も (32 ファイル)。
- `two_tier_balance=.true., two_tier_level=7` と `two_tier_level=6, two_tier_test_shift=1` (`-DTWO_TIER_CHECK` ビルド):
  完走、毎 step の balance で `TWO_TIER STATUS ... ok`、`check_particle_tree` 0 件、footprint 1148 区間。
  hydro の総和は単一分割と 1e-15 の相対差 (集団和の順序)。LIGHT ビルドでも同じ namelist で完走。
- restart: 2 層の output から `.true.` で再開 (Ls と footprint を読む)、`.false.` では
  "Restart file has a two-tier decomposition" で停止、単一分割の output から `.true.` で再開して最初の balance で 2 層に移行。
- 本番規模 (ks の 448 / 1792 ランク) での効果は ks の 7g の結果 (355 s/step vs 415 s、メモリ 20%) に拠る。こちらでは未実施。
