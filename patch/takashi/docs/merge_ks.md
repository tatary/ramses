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

## Ls の自動選択と Hilbert key の高速化 (2026-11-17): `two_tier_level=-1`, `patch/takashi/hilbert.f90`

`two_tier_level=0` の経験則 (細かい層の oct が全体の `two_tier_fine_frac` 以下) は ks でも使われていない (本番は固定 13)。
`two_tier_level<0` で、balance のたびに候補 Ls (levelmin 〜 grid のある最も細かい level) を全部切って採点する
(`load_balance.f90` の `two_tier_auto_ls`; 切り方は本番と同じ `two_tier_cut_sorted`、key のソートは 1 回)。

- 時間: critical path = Σ_level niter × (仕事量 (セル×80+粒子) のランク最大) を oct 数だけの理想値で割ったもの
  (ks の `lb_two_tier_estimate.py` と同じ考え方)。層の食い違い ghost (level Ls+1 の grid の所有者 ≠ 親 grid の所有者) は
  `two_tier_auto_ghost` (既定 0.5) 個分の所有 oct として親の所有者に数える (親は 1/8 個分を子の所有者に)。
- メモリ: ランク最大の max((所有 oct + 食い違い ghost)/ngridmax, 粒子/npartmax)。通常の表面 ghost は数えない。
- 選択: メモリが `two_tier_auto_memmax` (既定 0.7) 以下の候補のうち時間最良 (どれも超えるならメモリ最小)。現在の Ls は
  時間が最良の `two_tier_auto_tol` (既定 0.05) 以内でメモリに収まる限り保持する (Ls を変えると層が入れ替わる level の
  grid が全部移動する)。メモリを同点決勝に使う版 (「時間が 5% 以内のうちメモリ最小」) は、粒子数のランク最大が
  balance ごとに数ポイント揺れるせいで Ls が 11↔12 と振れたのでやめた。
- log に候補ごとの表 (`Two-tier auto Ls`) が出る。

### `hilbert3d` (QUADHILBERT) の key の組み立て

採点 (0.6 s/回) と cpu map (3.1 s/回) の内訳を取ると、ほぼ全部が `cmp_ordering` → `hilbert3d` だった (1 key 約 8.5 µs:
cpu map のセルのループ 1.95 s、footprint 0.45 s、item 0.24 s、採点の key 0.48 s; 切る bisection は 0.06 s、ALLREDUCE は
1 ms)。原因は key を実数にする最後のループが bit ごとに `real(2,kind=16)**i` (ソフトウェアの 4 倍精度べき乗) を呼ぶこと。
`patch/takashi/hilbert.f90` (上流のコピー) で、QUADHILBERT のとき bit を整数 2 個 (0-59、60-) に詰めて 1 回だけ 4 倍精度に
直す形にした (値は同じ: 96 bit までの整数は 4 倍精度で正確)。QUADHILBERT なしのコードは不変。
境界の展開/縮小 (`authorize_fine` の `cmp_minmaxorder`) も同じ所を通るので一緒に速くなる。

### 検証

MUSIC zoom、`kstest_long` の output_00003 (step 108、level 16 まで 33 万 oct) から 12 ランクで再開、nremap=1
(`claude_scratch/run/tt*`、バイナリ `claude_scratch/build_tt_new` / 変更前 `build_tt_base`)。

- bit 一致: 固定 Ls=13 の 4 step 後の 72 ファイルが変更前バイナリと一致 (hilbert の変更込み)。全体 221 s → 65 s。
- load balance 1 回 (Ls=8): cpu map 3.2 → 0.16 s、展開 1.0 → 0.07 s、縮小 0.8 → 0.08 s。採点 0.6 → 0.06 s。
  step は 15–17 s → 9.7 s、load balance の割合は 4–7 割 → 約 1 割。
- hilbert 修正後の固定 Ls との比較 (落ち着いたあとの step 時間):

| Ls | step (s) | 採点 (ghost 重み 0 / 0.5 / 1) |
|---|---|---|
| 8 | 9.7–9.9 | 1.94 / 2.10 / 2.08 |
| 自動 (既定) → 9 で安定 | 9.2–9.7 | 2.10 |
| 12 | 10.2–10.6 | 2.0 / 2.2 / 2.3 |
| 13 | 10.8–14 (6 step では落ち着かない) | 2.07 / 2.4 / 2.8 |
| 14 | 10.1–10.3 | 1.68 / 2.2–2.6 / 2.9 |
| 15 | 11.1–11.3 | 2.2 / 2.3–2.9 / 2.7 |
| 16 (細かい層なし) | 11.1–11.4 | 2.4 / 2.4–3.0 / 2.7 |

- 「最速」と「選ばれた Ls」が食い違って見えた理由 (hilbert 修正前: Ls=8 15.3 s、自動 Ls=11 17 s): ソルバの時間は同じ
  (10.5 s/step) で、差は全部 load balance (5.6 vs 7.0 s/回)。採点は load balance 自体のコストを見ておらず、当時は
  それが壁時計の半分だった。修正後は差が消えた (9.7 vs 9.8 s)。
- ghost の重み: 0 だと Ls=14 が最良に見えるが実測は Ls=8 より 4% 遅く、1 だと Ls=14 が最悪に見えるが実測は 12 より速い。
  間の 0.5 を既定にした。実測の幅が全体で 15% しかなく採点の分解能 (数 %) と同程度なので、これ以上は詰められない。
- DICE ビルド (本番設定、`dice_rt_amr` 4 ランク 6 step): 既定 (2 層なし) と `two_tier_level=7` が変更前バイナリと bit 一致
  (24 ファイルずつ)。`two_tier_level=-1` は完走、balance 11 回 `TWO_TIER STATUS ok`、Ls は初期の balance で 6→7→8 と動いたあと
  6 で安定 (時間スコアは 6–8 で 1.14–1.16 とほぼ同じ)、採点 0.011 s/回。
- 未確認: 本番規模 (448/1792 ランク、level 22; 時間スコアの差がもっと大きいはず) で選ばれる Ls と ks の固定 13 の比較。

## 星形成の移植 (2026-11-13): `sf_fg` と `star_formation_fg.f90`

### 方式

- ks の `star_formation.f90` (7500 行: 局所 Pop II、SFC、Pop III、非局所枯渇) は **別ファイル**
  `patch/takashi/star_formation_fg.f90` に `star_formation_fg(ilevel)` として持ち込み、`&SF_PARAMS` の
  `sf_fg=.true.` で `amr_step` が呼び分ける。takashi の `star_formation` (`sf_model`) は手つかずで、既定は bit 一致。
  ks の `contains` 構造 (ホスト変数を共有する 30 本余り) をそのまま保てるのでこの方が簡単だった。
- **粒子 family は `sf_fg` で切り分け**。`pm_commons` に `FAM_POP3=3, FAM_SN=4, FAM_DEAD=5, FAM_BH=6, FAM_SFC=7,
  FAM_PSC=8` (番号は ks と同じ → ks の part ファイルから restart 可) と `is_pop3` 等を追加。上流の
  `FAM_CLOUD/DEBRIS/OTHER` (3-5) と番号が重なるので、`sf_fg=.true.` では sink (cloud) と `f_w>0` を
  `read_hydro_params` で禁止。`NFAMILIES=8`。header の family 表は Mera が位置で読む 11 行の並びを維持し、
  `sf_fg` では 3-5 行目の名前を `pop3/supernova/dead_pop2` に、6-8 (`bh_pop3/SFC/PSC`) は fields の後に別ブロック。
  名前に空白を入れない (Mera/yt が `split` で読む)。
- **ログは takashi の仕組みに統一**: `sf_log.f90` (module `sf_log`: `sf_log_open/close/write_star`) を新設し、
  `star_formation`, `feedback`, `star_formation_fg` が共用。ランクごとの `stars_XXXXX.outNNNNN`、restart は
  `output_restart` のスナップショットで巻き戻す従来方式。列は従来 + 末尾に `family`。ks の
  `write_log_star_formation` (rank 1 への MPI 集約、`Nnew_max=1e5` の自動配列) は削除し、粒子生成の 5 か所
  (局所 Pop II、SFC 星、SFC 印、Pop III ×3) で直接書く。`logSFC` (SFC 1 個 1 行、rank 1) は ks のまま。
- `my_utils.f90` (ufront、位置→cpu/cell)、`sbs_imf.f90` (star-by-star IMF) は ks からほぼそのまま
  (`include 'mpif.h'` → `use mpi_mod`)。`ufront` (uold と同じ大きさ) は `init_amr` では確保せず
  `set_ufront` の初回に確保 (`ensure_ufront`) — 非局所枯渇か SFC を使うときだけメモリを食う。`flag3` は
  `sf_fg` のとき `init_amr` で確保。
- `output_analysis.f90` は **stub** (ugvtk 出力: `output_around_SFC/star[_at_given_age]`)。持ち込むかは未決。
- `cr_write_patch.sh`: `cat -e` が UTF-8 コメントを `M-"` に書き換えて埋め込み文字列を壊していたので `sed` に。

### ks からの変更点

- 全 level で呼ばれるが nlevelmax 以外は即 return (ks と同じ)。`amr_step`/`rt_step` の `uold(1:4)` 非同期交換の
  start は `sf_fg` では nlevelmax だけにし、FG 側で end して**仮想セルを速度に変換**する。ks は変換していないので
  領域境界の div v / 27 セル速度分散で仮想セルの運動量を速度として使っていた (ks 側のバグ; 本番は 448 ランク)。
- `boxlen /= 1` 対応: ks の球殻・球内操作 (`get_shell_quantities`, `modify_sphere_*`, `crop_and_adjust_SFC`,
  `distribute_stars`) はセル体積を `(1/2^l)^3` (単位箱) で計算し、境界判定を `[0,1]` としていた。
  `vol_cell=(dx_cell*scale)**3`、境界は `boxlen` に。cosmo (boxlen=1) では同じ。
- `mp0`/`msp0`: ks の `if(star_by_star.or.fire3_yields) mp0=mp` を `INIT_STELLAR_MASS` / `STELLAR_POPULATION_MASS`
  の ifdef に。takashi の feedback (msp0) がそのまま FG の Pop II に効く。
- 削除: test problems (`is_test`, `test_problem`, `add_star_particle_for_test`)、`fof_commons`、デバッグ用の
  自己遮蔽 print (`comp_SH2`, `getMu`, `cmp_LVG_length_scale`)、BH/SN 関連の設定 print (別段階)、
  `PXFMKDIR` → `execute_command_line`。
- 大きな固定長配列 (`Nrbin×Npop3_max/SFC_n_max/Nnew_max` の殻配列、最大 50 MB) は allocatable に
  (ifort ではスタックに載って落ちる)。
- 末尾の `make_virtual_fine_dp` × nvar × level は `make_virtual_fine_dp_nv` (level あたり 1 メッセージ)。
- **Pop III の雲半径で止まらない** (2026-11-19、`get_mdot_dependent_pop3_mass`): ks は (a) 流入している殻
  (VelR < 0) が 1 つも無い、(b) 選んだ雲半径が箱の境界に掛かる、のどちらでも `stop` した (ランク 1 だけが stop)。
  (b) は箱からはみ出す殻 (`get_shell_quantities` は周期の折り返しをしないので質量が欠けている) を探索から外し
  (`shells_inside_box`)、雲半径が無ければ警告を出して Mdot = 0 → `Mpop3 = 0` (その step は星もガス枯渇も無し;
  `pop3_model=6` の質量下限と同じ経路)。止まらずに走っていたランの結果は変わらない。殻の最大半径は
  1000 dx_min なので (b) は levelmax の低い箱 (DICE テスト: 箱の 4 倍) でだけ起きる。(a) は `loc_vir` の
  virial 条件を課していればまず起きない。
- **`loc_vir` の速度分散は 1 次元** (2026-11-21): ks は 27 セルの立方体の 3 次元の分散 (x, y, z の分散の和) をそのまま
  等温音速の 2 乗 (1 次元の量) に足していた。`sigma_1d^2 = sigma_3d^2 / 3` に直した (Pop III と Pop II / SFC 印の両方の判定)。
  同じ `loc_vir` でも条件が緩くなる (乱流が支配的なセルで alpha_vir が最大 1/3)。`loc_vir <= 0` のランは変わらない。
- `SFdata/pop3prof.*.dat` の Mdot の単位 (2026-11-19): ks はヘッダの `Mdot_cloud` と `Mdotshell` 列を `[Msun/yr]` と書きながら
  Myr で換算していた (値は Msun/Myr)。Msun/yr に直した (標準出力は元から Msun/yr)。ks の古い `pop3prof` と比べるときは 1e6 倍違う。
- **SFC 印 (`pop2_SFC_marker`) の `logSFC` / `logSFC_1Dprof`** (2026-11-19、`create_SFC_markers`): 印には SFC の
  半径・質量が無いので、Pop III の質量に使う雲 (流入している殻のうち Menc/M_BE が最大のもの; 殻の計算は
  `shell_derived_quantities` / `cloud_from_shells` で Pop III と共用) を書く。`logSFC` は列の位置を SFC と同じにして
  R_cloud, M_cloud, N_star = M_star = 0, 雲内の平均 nH と Z, pos、最後に **Mdot_cloud [Msun/yr]** の列を足す
  (ヘッダも変える; SFC モデルでは `pop2_SFC_marker` が切られるので 1 つのファイルに混ざらない)。
  `logSFC_1Dprof` は SFC と同じブロックの後に `rcen, VelRshell [km/s], MBEshell, Menc/MBE, Mdotshell [Msun/yr]` の行
  (セルの無い殻は -1、Mdot は 0; `Nshell` で区別する)。雲半径が無い印は `ir = -1` で R, M, Mdot が 0。
  コストは印ができた step に全 leaf セルの掃引 1 回 (`get_shell_quantities`)。

### 検証 (dice_rt_amr の 6 step、4 ランク; `tests/dice_rt_fg.nml`)

| テスト | 結果 |
|---|---|
| 既定 (`sf_fg=.false.`) | hydro/part/grav/rt/amr 46 ファイルが変更前と bit 一致。`stars_*.out` は family 列を除き同一 |
| 局所 Pop II (`pop2_model=2`) | 2675 セルで 2700 粒子 (100 Msun)。log の nH 前後から Σ 枯渇質量 = Σ 星質量 |
| 同、LIGHT_MPI_COMM ビルド | heavy と 46 ファイル bit 一致 |
| 同、`nonlocal_depletion=.true.` | 完走 (`ufront` 確保、球からの枯渇) |
| 同、restart (output_00001 から) | 完走、時刻一致。bit 一致ではないが既定の restart も同様 (既存の性質) |
| Pop III (`Zpop3=10`, `pop3_model=1`) | 10104 候補 → 120 Msun の Pop III、family 3 |
| Pop III `pop3_model=6` (mdot 依存) | 1D プロファイルと `SFdata/pop3prof.*` まで動き、雲半径が箱 (100 kpc) を超えて ks の境界チェックで停止 (粗い試験箱のため) |
| SFC (`pop2_model=4`) | 12 個の SFC (R 1-1.4 kpc、M 4-6e6 Msun、f* 0.045-0.05)、mass-weight 配置、crop、`logSFC` |
| SFC + 二層分割 (`-DTWO_TIER_CHECK`) | 完走、STATUS ok。SFC の個数は同じで星の割り当ては分割依存 (候補の除外順が cpu 順なのは ks の性質) |
| ks 設定 (QUADHILBERT, NGROUPS=5) | ビルド成功 |

### 残り (TODO 3b)

- ugvtk 出力 (`output_analysis`、3600 行) の移植可否。
- `SFC_flg_thin` の cooling 側 (`collect_active_SFC` による遮蔽解除)、`flg_thin_star_cell`。
- Pop III / SFC は takashi の feedback・RT では `is_star` でないので放射も SN も出ない (フィードバック段階で ks の
  `rt_spectra`/`feedback` を対応させる)。
- ks の part ファイルからの restart: ks は `mp0` を `star_by_star.or.use_initial_mass` のときだけ書き、takashi は
  `INIT_STELLAR_MASS`/`STELLAR_POPULATION_MASS` で `mp0`,`msp0` を読む。ks 出力を読むビルドでは
  `STELLAR_POPULATION_MASS=0` にするか reader を切り替える必要がある。
- stdout が星 1 個ごとに数行 (ks のまま)。本番で抑える手段 (`verbose` ゲート) を検討。
- `random_value_0to1` (SFC の配置、IMF 抽出) の乱数状態は restart で保存されない (ks と同じ)。

### 非局所枯渇の一括処理 (2026-11-13、`nonlocal_depletion_batch`)

ks の `deplete_gas_from_spheres` は新星 1 個ごとに「全 level の全セル走査 + ALLREDUSE 数本」を逐次に行う
(重なった球で引きすぎないため、前の星の枯渇を反映してから次の星の殻質量を計算する)。1 step に N 個できると
集団通信が O(N) 回になるので、`deplete_gas_from_spheres_batch` (既定 `nonlocal_depletion_batch=.true.`) に置き換えた。

- 新星リスト (位置・質量) を 1 回共有 (`collect_newstar_info`、ks のまま)。
- 各星について octree を根から辿り (`sphere_walk`/`sphere_walk_grid`)、球と交わる立方体だけ降りて、自分の
  active grid (`cpu_map(father(igrid))==myid`) の leaf セルを集める。殻質量は最初 `rbin(13)` (≈4 dx_min) まで、
  足りない星だけ範囲を倍にして再計算。ALLREDUCE は `(Nrbin × N)` 配列 1 本/round。
- 半径と `f_i = M_i/M_enc,i` は ks と同じ規則 (`M_enc > M_i/0.2` の最小 rbin)。
- 適用は**加法**: 球に含まれるセルは `ρ → ρ (1 − Σ_i f_i)`。引く質量は厳密に `Σ M_i`。重なったセルで
  `Σ f_i > f_sum_max (=0.5)` になったら、その球の星の半径を 1 bin 広げて (f が減る) 検査し直す (反復ごとに
  N 整数の ALLREDUCE 1 本)。`Nrbin` まで広げても超える場合だけ `f_sum_max` で clamp し、引けなかった質量を
  警告に出す。重なりがなければ ks の逐次版と同じ結果 (fg4 で hydro の総和が一致)。
- `uold/unew/ufront` の更新は `modify_sphere_by_factor` と同じ (nlevelmax は密度のみ、粗い level は全変数と unew の差分)。
- 試験: `f_sum_max` を一時的に 0.05 にしたビルドで反復経路 (2–5 回、clamp なし、`M_removed = M_star`) を確認。
  heavy/LIGHT bit 一致、二層分割 (`-DTWO_TIER_CHECK`) で完走。
- `nH_postSF>0` の版と Pop III の `deplete_gas_around_pos` は逐次のまま。

## フィードバックの移植 (2026-11-13)

ks 本番 (`ngridmax_fire3.nml`) は `thermal_feedback` + FIRE-3 yields + **SN reservoir** + 冷却停止タイマーで、
`mechanical_feedback` / `spread_feedback` / `nonlocal_SN` は使っていない (ks の amr_step にも「このプロジェクトでは使わない」
とある)。takashi の `feedbk` は既に FIRE-3 (msp0 方式) なので、差分だけを takashi の `feedback.f90` に足した。

### 入れたもの

- **Pop III の SN** (`feedbk`、`sf_fg` の `FAM_POP3` 粒子): Schaerer (2002) の寿命 (`cal_tpop3`; pop3_model 1/2 は 4 Myr) の
  経過時に 1 回だけ爆発。エネルギーは Wise+2012 の質量依存式、pop3_model 1/2 は `pop3_ESN` (本番 3e52 erg)。残骸質量・金属
  は ks の表 (`cal_mrem`, `cal_mmet`; 1/2 は残骸 100 Msun・金属 9 Msun)。300 Msun 超 (model 2 の 1000 Msun) は直接崩壊。
  SN 後の family は `FAM_BH` (sink の段階で bhpop3 の種)。ログは `stars_*.out` の event=1 行 (family 3 で記録してから BH に)。
  `thermal_feedback` の粒子選択は `is_star .or. (sf_fg .and. is_pop3)`。
- **SFC → PSC** (`thermal_feedback`、`sf_fg`): `SFC_t_life` を過ぎた `FAM_SFC` を `FAM_PSC` に。
- **SN reservoir** (`sn_sphere_stats.f90`, `sn_feedback_radius.f90`, `sn_energy_reservoir.f90`; ks からほぼそのまま、`use mpi_mod`):
  `&FEEDBACK_PARAMS` の `SN_reservoir=.true.` で、feedbk は純 SN 熱 (`f_esn × E_SN × N_SN`、風の熱は除く) を unew に入れず
  `pending_heat` として除き、`thermal_feedback` 末尾で全ランク共通のイベント列 (位置・E・idp、ALLREDUCE 3 本) を作って
  半径探索 (冷却半径 `28.4 pc nH^-3/7 (E/1e51)^2/7 f(Z)^-1/7` と `r < r_cool/2`、候補は host 半対角から 1.2 倍刻みで
  `Ncell_SN_max × dx` まで、集団通信 2 本) → `deposit_sn_reservoir` が球内の自分の leaf に体積比例で `sn_pending` [erg/コード体積]
  を積む。`amr_step` の refine 前の `prepare_sn_reservoir(ilevel)` が nlevelmax→ilevel の順に `uold` へ反映 (`upload_fine` と
  ghost 更新込み)。refine/derefine (`refine_sn_reservoir`/`coarsen_sn_reservoir`)、load balance の ghost 共有、defrag の並べ替え、
  hydro checkpoint 末尾の版付きトレーラ (reservoir 使用時のみ書く; 無いときは旧 checkpoint 扱い) を takashi の各ファイルにフック。
  粗い step ごとに rank 1 が `SN reservoir: released X erg into N cells` を出す (追加)。
- **冷却停止タイマー** `fcool_stop` (`&FEEDBACK_PARAMS`): `ivar_coolstop=nvar` の passive scalar に `fcool_stop × dx/c_s` [Myr] を
  置き、`cooling_fine` が正の間は冷却せず dt を引く。reservoir 使用時は反映時に host セルだけ、即時注入時は feedbk が設定。
  NPSCAL を 1 つ増やしたビルドが必要 (ks 設定は NPSCAL=3; `read_hydro_params` が検査)。
- **Pop III の放射** (`rt_spectra.f90`): `is_pop3` を線源に含め、Schaerer (2002) 表のフィット (`get_pop3_Nphoton[_pop3_model]`、
  群のエネルギーで帯を判定)、1e5 K 黒体の群特性 (`get_pop3_props`)、寿命 (`get_pop3_ageGyr`; model 1/2 は SN まで)。
  `&RT_PARAMS` に `p3_noFUV`, `p3_noEUV`。X 線/hypernova 群 (`isXray`, `rt_hypernova`) は未移植。
- `cooling_fine.f90`: `cooling_on` の宣言が `#ifdef RT` の中にあって RT=0 ビルドが落ちていた (cooling 移植時の退行) のを修正。

### 持ち込まなかったもの

`mechanical_fine[_pop3]`, `spread_fine`, `feedback_flag`, `nonlocal_feedback` (nonlocal_SN, nHmax_preSN)、star-by-star の SN
(Inoue+2011 残骸、`FAM_DEAD` 化)、`pop2_SN_model` (fixed/stochastic/mass_dep; takashi は常に FIRE-3)、`logSN` (takashi の
`stars_*.out` event=1 行で代替)、hypernova X 線。

### 検証 (dice_rt_amr 14 step、4 ランク)

| テスト | 結果 |
|---|---|
| 既定 (`sf_fg=.false.`, reservoir off) | 46 ファイル bit 一致 (hydro のトレーラは reservoir 使用時のみ) |
| Pop III SN 即時注入 (Zpop3=10, pop3_model=1, 1e50 erg) | 846 個が step 13 (4 Myr) に爆発、BH(PopIII)=846、event=1 行 family 3 |
| 同 + reservoir (3e52 erg) | 半径探索 846 件 (球 69、host 未解像 104、候補なし 673)、放出エネルギー = ΣE_SN (2.538e55 erg) |
| 同 + `fcool_stop=8` (NPSCAL=2 ビルド) | `ivar_coolstop=11`、完走 |
| Pop II FIRE-3 + reservoir | 1916 SN (全て host 未解像 = 単一セル)、放出 = ΣE_SN、event=1 行 family 2 |
| LIGHT / 二層分割 (`-DTWO_TIER_CHECK`) | heavy と 69 ファイル bit 一致 / 完走、STATUS ok、同じ放出量 |
| restart (output_00002 から) | 同じ 846 SN と放出量 (トレーラ読み込み経路) |
| RT=0 ビルド (sedov) | ビルド・実行 OK |
| ks 設定 (QUADHILBERT, NPSCAL=3) | ビルド OK |

## sink/BH の移植 (2026-11-13): `sink_particle_fg.f90`

### 方式

- ks の sink は上流より古い系統 (cloud 粒子を `idp=-isink` で識別、family なし、`mseed/acc_rate/acc_lum/dt_acc/level_sink`、
  Pop III 残骸 (`FAM_BH`) から sink を作る `make_sink_from_pop3_new`、BH_model 0-4 の降着・光度、ks 形式の `sink_XXXXX.out`)
  で、takashi が使う上流の `pm/sink_particle.f90` (clump ベースの生成、`FAM_CLOUD`、msmbh、stellar object) と同名ルーチンが
  中身違いで並ぶ。星形成と同じく **別ファイル + `sf_fg` で振り分け**: ks の `sink_particle.f90`/`init_sink.f90`/`output_sink.f90`
  (+ `load_balance.f90` の `check_cloud_count`) を scratchpad の `make_sink_fg.py` で変換して `sink_particle_fg.f90` (約 4000 行)
  にした。全ルーチンに `_fg` を付け (`create_sink_fg`, `init_sink_fg`, `backup_sink_fg`, ...)、`mpif.h` → `use mpi_mod`、
  `is_test`/clump/sgagn/`weightp`/`p3idpsink`/`restart_files` を削除。上流版はそのままリンクされ `sf_fg=.false.` で使われる。
- 振り分け箇所: `amr_step` (update_cloud, collect_acczone_avg, update_sink ×2, grow_sink, create_sink と出力後の
  `acc_rate=0`)、`init_part` (`init_sink_fg`)、`particle_tree` (`kill_entire_cloud_fg`/`create_cloud_from_sink_fg`)、
  `output_amr` (ks と同じく全ランクが同じ内容の ks 形式バイナリ `sink_XXXXX.out<rank>` / `sink_restart.out<rank>` を書き
  (restart では各ランクが自分のを読む)、rank 1 が yt 用の `sink_file_descriptor.txt` と ks のテキスト表 `sink_XXXXX.info`
  (`output_sink_fg`)、上流形式の csv を書く。一度 rank 1 だけが書く形にしたが、yt の RAMSES フロントエンド
  (`SinkParticleFileHandler`) は `setup_handler` で全ドメインのファイルの存在を要求する (中身はドメイン 1 しか読まない) ので
  戻した (2026-11-14、base 4.3.1 / main / `/mnt/synology2/takashi/smbh/ana/yt` の dev 版とも同じ)。yt には上流 csv を読む
  `SinkParticleFileHandlerCsv` (ptype `sink_csv`) もある)、`read_params`
  (`read_sink_params_fg`)、新たに patch に置いた `newdt_fine.f90` (降着 dt: `c_acc × min(M_gas, M_sink)/Mdot`、最高 level で
  `compute_accretion_rate_fg`) と `force_fine.f90` (`f_gas_sink_fg`: 共動の softening)、`to_dice/rho_fine.f90`
  (`no_sinkcloud_refine`)。
- **cloud 粒子の family**: sf_fg では 3 が Pop III なので、FG の cloud には `FAM_CLOUD_FG=9` (NFAMILIES=9、header の FG
  ブロックに `sink_cloud`) を与え、`is_cloud` が `sf_fg` を見て 3/9 を切り替える (pure 関数からモジュール変数の参照は可)。
  共有カーネル (`move_fine` の vsnew 上書き、`synchro_fine` の fsink 集計、`rho_fine` の `sink_refine`) はこれで動く。
  ks 側のルーチンは `idp<0` のまま (両方成り立つ)。
- `&SINK_PARAMS` (`read_sink_params_fg`): 上流の項目のうちこの `pm_parameters` にあるもの + FG の `bhpop3`, `BH_model`,
  `sink_mass_min_pop3BH`, `sink_mass_min_pop2BH`, `no_sinkcloud_refine`, `sgagn_T2max` (BH_model=3 の熱注入上限),
  `sink_check_vars` (`check_sink_vars_among_cpus_fg`; 既定 off)。`amr_parameters` に宣言。`nol_accretion`/`sink_drag`/
  `feedback_scheme`/`v_AGN` は落とした (旧 accrete_sink 専用)。
- BH の放射 (`rt_spectra`): `set_BH_SED` (α=-1.5 のべき、11.2 eV-10 keV; この版の積分器は帯が重ならないと止まるので
  `L_to_Ndot_BH=0` の群は飛ばす)、`update_SED_group_props` に acc_lum で重み付けた群特性の寄与、`star_RT_feedback` 末尾の
  `BH_RT_feedback` (sink の位置の leaf セルが ilevel のとき光子を入れる; 担当は `cmp_cpumap_lev` で求めた ilevel grid の
  所有者 = my_utils の `get_cpu_of_pos` と同じ)。`rt_star` の起動条件に `nsink>0 .and. BH_model 1/2` を追加。
- `init_sink_fg` は ks の配列に加えて共有コードが触る `msmbh/dmfsink/graddescent_over_dt` (0) も確保し、全部 0 初期化。
  `compute_ncloud_sink` は上流のものを呼ぶ。`get_BHL_accretion_rate_fg` は非 cosmo でも G=1 で動く (ks は停止)。
  `d_sink` が決まらなくても止めない (clump 生成でしか使わない)。

### 検証 (4 ランク)

sink 試験は `tests/dice_rt_fg.nml` の Pop III SN 設定 (Zpop3=10, pop3_model=1, 3e52 erg, 18 step) に `sink=.true.` と
`&SINK_PARAMS` (bhpop3, BH_model=2, sink_mass_min_pop3BH=10 Msun, bondi, ir_cloud=4, smbh) を足したもの。

| テスト | 結果 |
|---|---|
| 既定 (`sf_fg=.false.`, dice_rt_amr 6 step) | 変更前バイナリ (worktree) と 44 ファイル bit 一致 |
| FG 星形成のみ (dice_rt_fg, sink なし) | 同 bit 一致 (header に `sink_cloud 0` の行が増えるだけ) |
| sink 試験 (heavy) | step 13 の 846 SN → 粗い step 6 で 846 sink、合体 833 → 13 → 11; step 15 の 61 SN も sink 化して合体; cloud 粒子 11×2109=23199、cloud 数チェック 0 件不一致; Bondi 降着で 100 Msun 種が 1e4 Msun 級に (level 8 の粗い試験箱なので物理的な意味はない); `logsink<id>` 13 本、`sink_XXXXX.out<rank>`/csv/info |
| LIGHT_MPI_COMM | heavy と 44 ファイル bit 一致 (sink ファイル含む) |
| 二層分割 (`-DTWO_TIER_CHECK`, two_tier_level=7) | 完走、STATUS ok、cloud 数チェック ok; 総和は 1e-6 程度の丸め差、sink の質量・位置は合体順の違いで 2 個が僅差 |
| restart (output_00003, `sink_00003.out<rank>` から) | 完走、sink 11 個の質量・位置は 1e-9 で一致、星形成/SN イベント (stars ログ) は同一; hydro は丸め差 (FG/フィードバック段階の restart と同じで bit 一致ではない; 乱数状態と `acc_lum` などの派生量は保存されない) |
| ks 設定ビルド (QUADHILBERT, NPSCAL=3, DICE なし) | ビルド OK (`patch/takashi/rho_fine.f90` 経由) |

## ugvtk 出力の移植 (2026-11-14): `output_analysis.f90`

ks の `output_analysis.f90` (3600 行) のうち本番で使う部分を、生成スクリプト (scratchpad `make_output_analysis.py`、
ks のブロックを切り出して変換) で `patch/takashi/output_analysis.f90` (2300 行、以前は stub) にした。

### 方式

- 出すもの: `ugvtk/<prefix>.<nstep>.<level>.vtk` (legacy VTK binary、STRUCTURED_POINTS、各 level で `NI_ugvtk`^3 セル、
  `l_range` level 分の入れ子、rank 1 が `.tmp` に書いて rename)。フィールドは nH, T, vx/vy/vz, Z, yHI/yH2/yHII/yHeI/yHeII/yHeIII,
  Np_bin1..nGroups (RT)、phi_grav、rhoDM/rhostar。中身は `output_uniform_data` (ks のまま: `set_ufront(lev,2)` と
  無条件の `set_ufront(lev,3)`、`set_particle_density` の粒子密度、level ごとの ALLREDUCE で rank 1 に集約)。
- 経路: (1) SFC 形成時 `output_around_SFC` (`star_formation_fg`、`SFC_output_minmass`)、(2) Pop III 形成時 `output_around_star`
  (`p3_output_minmass`)、(3) `idp_ugvtk` の星の形成時と年齢ごと (1 kyr … 300 Myr の 33 点、`output_around_star_at_given_age`)、
  (4) 全 SFC の 0.1/1/3/5/10 Myr (`SFC_output_at_given_age`、重い)。(3)(4) は `amr_step` から毎 level 呼ぶ
  `output_analysis_data` が最も細かい level でだけ動かす (ks と同じ位置: levelmin の出力ブロックの直後)。
- **主スイッチ `ugvtk` (`&OUTPUT_PARAMS`、既定 .false.、TO ADDED)**: ks には無く常に有効だったが、既定だと Pop III 1 個ごとに
  数十 MB 書くので takashi では明示的に on にする。`idp_ugvtk`, `NI_ugvtk`, `SFC_output_at_given_age`, `p3_output_minmass`,
  `SFC_output_minmass` も `&OUTPUT_PARAMS` で読む (後の 2 つは `&SF_PARAMS` にもある; 同じ変数を 2 つの namelist 群に入れるのは可)。
- `phi_sync_for_output` (ks: poisson_commons + interpol_phi + amr_step + refine_utils): 各 level の最初の subcycle で phi を保存し、
  出力の重力ポテンシャルを level 間で同じ時刻に揃える。この module の配列にして `ugvtk` のときだけ確保、`amr_step`
  (フラグ `flg_phi_sync_for_output` と `save_phi_sync_for_output`) と `refine_utils` (新 grid は父の値、消す grid は 0) にフック。
- `rhoprt` (ks: amr_commons、常に ncell×2) も module の配列で初回使用時に確保。`nstep_restart` を `amr_commons` に足し
  `init_amr` の restart で設定 (restart 後、捨てた枝の nstep 以降の vtk を `ugvtk/orphan_<nstep>/` に退避する ks の処理)。
- Pop III 形成時の出力 (pop3_model 1-3) と `idp_ugvtk` の星の形成時出力は、ks では `write_log_star_formation` の末尾
  (全ランクに集めた新星リストのループ) でやっていた。takashi は星ログを生成箇所で直接書く方式にしたのでリストが無い。
  代わりに `ugvtk_note_new_star(ipart)` (生成箇所 8 か所、rank 内で Pop III と idp_ugvtk の星だけ控える) と
  `ugvtk_flush_new_stars(ufront_updated, pop3_output)` (ks が write_log を呼んでいた 4 か所の集団点: 局所生成の後、
  `add_pop3_particles` の後、SFC 印と SFC 星の後; ALLREDUCE で集めて全ランクが `output_around_star`) に分けた。
- Makefile: `output_analysis.o` は my_utils/poisson/hydro/rt の module を使うので RT の module 群の後 (sn_* の後) に移した。
- **単位**: ks は boxlen=1 (code 単位 = user 単位) 前提で、星の位置 (user 単位、DICE では ~50 kpc) を 0-1 の箱に当てはめて
  `Ni` が負になった。領域と標本点を user 単位 ([0,boxlen)) に統一し (`get_cell_index`/`cmp_cpumap` は user 単位)、level 幅だけ
  `scale_u=boxlen/nx_loc` で換算 (`output_around_SFC/star`, `output_uniform_data` の `bscale`)。ORIGIN/SPACING は user 単位 × scale_l で pc。
- **`my_utils` の `get_cpu_of_pos` のバグ修正** (同じ原因): 点を grid 中心に丸めるときの grid 幅が code 単位だったので、boxlen≠1 では
  点から 1/boxlen 以内の別 grid に落ちることがあり (cmp_cpumap は user 単位)、`get_active_cell_index_*` (SFC/Pop III 4-6 の配置、
  ugvtk の標本) が所有者を誤る可能性があった。boxlen=1 (ks) では変化なし。pop2_model=2 / pop3_model=1 の試験は bit 一致。

### ks からの変更点

- `ugvtk` スイッチ (上)。`fout_*` 系 (x 軸の線データ、中心箱、sink/halo 周り、密度最大・ポテンシャル最小) と `output_line_data`
  等は未移植 (ks 本番で全部 -1)。namelist に `fout_*` が残っていると読めないので `namelist/ks/ngridmax_fire3.nml` ではコメントのまま。
- `rt_c_fraction` はこの版ではスカラー (ks は level ごとの VSLA)。非 cosmo では `texp=t`, `texp_ini=0` (ks は Friedmann 表を無条件に引く)。
  RT なしビルドでは放射場を出さず、電離度は中性 (`#ifdef RT`)。
- Pop III 形成時の出力に `p3_output_minmass` を pop3_model 1-3 にも適用する (ks は 4-6 だけで、1-3 は `write_log` 経路で全 Pop III を出していた。
  ks 本番の `p3_output_minmass=1d4` と pop3_model=1 の 120 Msun だと、ks は全部出し takashi は出さない)。
- ks の `output_around_SFC` は star_formation 内の呼び出しを前提に nlevelmax の uold を原始変数とみなして保存量に戻して出力し、
  終わったら戻す。`output_around_SFC_at_given_age` (amr_step から、uold は保存量) もこれを通るので、出力の nlevelmax セルは
  `set_ufront(lev,2)` 経由なら誤変換になるが、無条件の `set_ufront(lev,3)` (unew) が後で上書きするので実害は出力には無く、
  uold には往復の丸め差だけ残る。→ 2026-11-21 に修正 (下の「ugvtk: 読む状態の修正と高速化」)。
- 年齢ごとの出力は 1 step に最初に該当した年齢 1 つだけ (`exit`) なので、fine step (試験では ~0.3 Myr) より短い 10 kyr, 100 kyr は飛ぶ (ks と同じ)。

### ugvtk: 読む状態の修正と高速化 (2026-11-21)

**読む状態 (ks 2f3541a の移植 + 残っていた 2 経路)**。ks は uold ベースの ufront (`set_ufront(lev,2)`) を必ず unew ベース
(`set_ufront(lev,3)`、「WARNING」で囲まれたループ) で上書きしていた。unew が状態として有効なのは level の流体更新から次の
load balance / refine までなので、`amr_step` の冒頭 (`set_unew` の前) に走る年齢出力は壊れた値を読む。takashi では ks より
ひどく、DICE の試験 (`idp425_*yr`) は 33 本全部が壊れていた (全セルで v = 2900–3700 km/s、T < 0、nH が 1 桁小さい; 移植時の
検証は NaN の有無しか見ていなかった)。
- 修正: module 変数 `ugvtk_state_uold`。`output_analysis_data` (amr_step) の中だけ真で uold (option 2、上位 level からの upload
  付き) を読み、star_formation の中からの呼び出し (形成時の入れ子) は従来どおり unew (option 3) を読む。`output_around_SFC` の
  uold の「原始変数 → 保存量」変換は `ugvtk_state_uold` が真のときはしない (amr_step では uold は最初から保存量)。
- ks 2f3541a は `use_current_ufront` の真偽で分けたが、それだと (1) SFC の年齢出力が `output_around_SFC` の変換で最細 level の
  速度・温度・金属・電離度を壊し、(2) star_formation の中の `use_current_ufront=.false.` の呼び出し (takashi では
  `ugvtk_flush_new_stars`) が原始変数の uold を読む。ks 側にも同じ修正を入れた (`Ramses-FG/patch/ks/NOTES_pop2_model_comparison_2026-11-11.md` 末尾)。
- 結果: 形成時の出力は変更前とバイト一致。年齢出力は変わる (DICE 試験: |v| 数十 km/s、T 1e4–6e5 K で、同じ step の形成時出力と
  nH が一致)。年齢出力のある step でも uold に丸め差が入らなくなった。

**高速化** (出力の中身は変えない; 形成時の出力で変更前とバイト一致を確認):
1. ufront の準備 (全 level の `set_ufront`、通信あり) を 1 ファイルごと → 1 イベント (星 / SFC 1 個の入れ子) に 1 回。さらに
   同じ瞬間に続けて出すループ (`ugvtk_flush_new_stars`、pop3_model 4–6 の形成時、SFC の年齢出力) は `ugvtk_batch_begin/end` で
   囲み、ufront と粒子密度 (`set_particle_density`) を全イベントで 1 回にした。SFC の形成時の 2 か所はループの中で粒子を
   作るので 1 個ごとのまま。
2. 点 → セルの対応: ks は全ランクが level ごとに全点について「担当ランクの計算 (Hilbert キー)」と「木の探索」をし、level ごとに
   立方体全体の `MPI_ALLREDUCE` を 4 回していた。今は 1 点につき木を 1 回たどり (ilevel まで)、着いたセルの格子の担当
   (`cpu_map(father(igrid))`) が自分なら自分の点とする。通信なし。
3. データの収集: 全ランクが立方体全体 (nvar + 光子 + phi + 粒子密度) を `MPI_ALLREDUCE` で足していたのを、自分の点だけを rank 1 に
   `MPI_GATHERV` で送る形にした。大きな配列は rank 1 だけが確保する。rank 1 は「どのランクも名乗らなかった点 / 複数が名乗った点」を
   数え、あれば WARNING を出す (試験では 0)。
4. ファイルの rename を `system('mv -f ...')` (1 ファイルごとに fork) から C ライブラリの `rename` に。
- rank 1 が 1 イベント (バッチ) ごとに ` ugvtk <名前>: N files; prepare / search / gather / write = ... s` を標準出力に出す。
- 手元の時間: DICE 4 ランク、Pop III 846 個の形成時出力 (2538 ファイル、`NI_ugvtk=4`) で star formation の行が 73.9 → 1.2 s
  (ラン全体 152 → 81 s)。星 1 個の年齢出力 (`NI_ugvtk=32`、36 ファイル) は 1.90 → 1.4 s (残りは 1 回分の準備)。ks zoom
  (12 ランク、`NI_ugvtk=64`、11 ファイル) は 1 イベント 1.7 s (準備 0.30 / 探索 0.27 / 通信 0.26 / 書き出し 0.83; 変更前は未計測)。
  ランク数に比例する無駄を削ったので本番規模で効くはずだが、クラスターでは未計測。
- 検証 (変更前のバイナリとの比較、vtk + hydro/part/grav/rt/amr): DICE Pop III 846 個 2578 ファイル一致、同 `two_tier_balance=.true.`
  2578 一致 (この設定では細かい側の格子が 0 で 2 段分割の検証としては弱い)、DICE 星 1 個 `NI_ugvtk=16/32` は高速化だけの版で 116 一致
  (修正後は年齢出力 33 本だけ変わる)、ks zoom 71 一致。

### ugvtk: ファイル名の番号を粒子の ID に (2026-11-21)

- 変更前 (ks と同じ): `SFC<n>` / `Pop3_<n>` の `<n>` は「その step に出力する SFC / Pop III の中での通し番号」。step ごとに 1 から
  振り直すので、形成時の `SFC1.<nstep>` と後の `SFC1_100kyr.<nstep>` が同じ SFC とは限らず、step 数から逆算するしかなかった。
- 変更後: `<n>` を粒子の `idp` にした。`SFC<idp>.<nstep>.<level>.vtk` (形成時; SFC と `pop2_SFC_marker` の印)、
  `SFC<idp>_100kyr.<nstep>.<level>.vtk` (年齢出力)、`Pop3_<idp>.<nstep>.<level>.vtk`。`<nstep>` (fine step の通し番号) は残す
  (restart 時の `orphan_<nstep>/` 退避が使う)。`idp<idp_ugvtk>…` は元から ID。桁は `I0` (幅制限なし)。
- 実装: SFC / 印の形成時は出力の直後に `nstar_tot+1` で粒子を作るので `output_around_SFC(nstar_tot+1, …)`。年齢出力は位置と一緒に
  `idp` を集める。Pop III は pop3_model 1–3 が `ugvtk_flush_new_stars` の控え (`note_idp`)、4–6 は `add_pop3_particles` の前に出すので
  同じ数え方 (`pop3_n_multi`: 単独 1、連星 2、多重 n) で先に ID を求める。連星・多重星は最初 (最も重い) の星の ID。
  星ができなかった候補 (質量 0: 流入殻が箱の外、多重星の最小質量未満) は ID が無いので出さない (変更前は
  `p3_output_minmass<0` のときだけ `Pop3_<n>` として出ていた)。
- ks 側 (Ramses-FG) は変えていない (FG 系のランはこちらのコードでやり直す予定)。ks の出力と比べるときは名前が合わない。
- 検証 (DICE 4 ランク、変更前のバイナリと比較): Pop III 846 個 (pop3_model 1) は 2538 本の中身が同じで ID は `New star added` の
  ログと一致; pop3_model 6 (候補 90、うち星ができた 84 系、2–5 重) はファイルの ID = 各系の 1 番目の星の ID で、84 系 252 本の中身は変更前の対応するファイルと同じ; SFC (`pop2_model=4`、
  `SFC_output_at_given_age`、8 step) は `SFC1 → SFC1_100kyr → SFC1_1Myr → …` と同じ ID で追える; ks zoom (12 ランク) は 11 本同じ。

### 検証 (4 ランク、sink 試験の namelist `tests/dice_rt_fg.nml` 相当、levelmin 6 / levelmax 8)

| 試験 | 結果 |
|---|---|
| `ugvtk=.false.` (既定、18 step) | 完走、変更前 (fs2) と 112 ファイル bit 一致 |
| `ugvtk=.true., idp_ugvtk=425, NI_ugvtk=16, p3_output_minmass=1d9` (18 step) | 完走 (186 s、off と同じ)、`ugvtk/idp425.<nstep>.{6,7,8}.vtk` (形成時) + `idp425_1kyr/_1Myr/…/_10Myr` 各 3 本 = 36 本 13 MB; 寸法 16³、原点 = 星の位置 − 8Δx、22 フィールド NaN なし (HII 領域で yHII→1、光子密度は星の周りで非零、phi_grav は phi_sync 経由、rhostar=0 は FAM_STAR 無し)。hydro/part/grav/rt の出力は off と bit 一致 (出力が物理を変えない) |
| `ugvtk=.true., NI_ugvtk=4, p3_output_minmass=-1` (7 step、Pop III 846 個の形成時出力) | 完走 (158 s)、`Pop3_<k>` (現在は `Pop3_<idp>`) 846 個 × 3 level = 2538 本 20 MB、`.tmp` 残りなし |
| RT=0 ビルド | コンパイル OK (放射場なし、電離度は中性) |

## ks との比較ラン (2026-11-14): `namelist/ks/l5_14_test.nml`

ks 版バイナリと移植後のバイナリを同じ cosmological 初期条件で走らせて比べた (方針 5 の「FG 初期条件での比較ラン」)。

### 設定

- 初期条件: ks の two-tier 試験に使った MUSIC zoom `work7/ramses_test/data/music/haloD_l5_14_z9` (level 5-14、Ωm=0.27、h=0.7)。
  ks 側の run は `work7/ramses_test/runs/base_np16` (ks `jobfiles/ws/two_tier/run_l5.sh` が `namelist/first_galaxy/l5_14.nml` に
  当てたもの: levelmin 5 / levelmax 16、16 ランク、6 coarse step、foutput=2、pop2_model=5 (SFC)、pop3_model=6、sink bhpop3)。
- ks 側: `Ramses-FG/bin/ramses0X3d` (dev_to、2026-11-11、-O2 -fcheck=bounds)、`runs/base_np16/namelist.nml` に
  `flg_voffset_ini=.false.` (初速オフセットは未移植) と `fof_find=.false.` (この設定では物理に効かない) を当てたもの。
  実行: `work7/ramses_test/claude_scratch/run/ksref`。
- takashi 側: `build_ks` (NGROUPS=5, NPSCAL=3, QUADHILBERT, DICE 空, -O3) + `namelist/ks/l5_14_test.nml`
  (`ngridmax_fire3.nml` と同じ翻訳規則で `base_np16/namelist.nml` を翻訳)。実行: `.../run/kstest`。
- 比較: 各 coarse step のエネルギー (log)、level ごとの grid 数 (log)、`tests/hydro_stats.py` (全セルの総和/min/max)、
  粒子数 (header)。スクリプトは scratchpad `cmp_ksreg.py`。

### 合わせるのに必要だった 2 点 (どちらも ks の挙動を移植)

1. **`omega_b`**: ks は `&RUN_PARAMS` の `omega_b0` (既定 0.0493) を `read_params` で無条件に `omega_b` に代入する。
   takashi は上流の `omega_b` (`&INIT_PARAMS`、既定 0.045)。試験 (0.044) では密度が一律 2.3% ずれていた。
   namelist で `omega_b` を合わせる (本番の翻訳 `ngridmax_fire3.nml` にも `omega_b=0.0493` を追加)。
2. **初期電離度**: ks の `init_flow_fine` は grafic に電離度ファイルが無いとき、cosmo なら Peebles の残留電離
   x_e = 1.2e-5 √Ωm/(h Ωb) (= 2.02e-4) を RT の電離度の既定値にする (上流は ATON のときだけ)。さらに ks は `rt_is_init_xion` の
   既定が .false. (上流は .true. で、`init_refine` が T_start の衝突平衡 x_e ~ 1e-20 で上書きする)。
   `patch/takashi/init_flow_fine.f90` (hydro からの新しいコピー) と `to_dice/init_flow_fine.f90` に同じ既定値を入れ、
   namelist で `rt_is_init_xion=.false.` を指定する (両 ks namelist に追加)。H⁻ 経路の H2 生成は残留電子で決まるので
   初代銀河の計算では本質的。非 cosmo (DICE) には影響しない (fs3 の回帰で確認)。

### 結果 (6 step、16 ランク)

| 量 | ks vs takashi |
|---|---|
| 初期状態 (output_00001) の密度・運動量・xHI・xHII | 全桁一致 (1.79M DM 粒子も同数) |
| 初期温度 T/μ | 128.58 vs 127.90 K (0.5%; `rt_evol_single_cell` の冷却率テーブル・定数の差) → 内部エネルギーの差 0.5% が以後も残る |
| grid 数 (level 1-14) | 全 step で一致 (level 14 が step 2, 5 で 5 個・2 個違うだけ、22 万個中) |
| epot / ekin / eint (step 1-6) | takashi の log の 3 桁の範囲で一致 |
| 密度の総和 / 運動量の総和 (output 2-4) | 1e-6 / 3e-5 の相対差 |
| xHII | 3e-4 の相対差; xHeII/xHeIII は下限値の差 (ks `x_MIN=1e-12`、takashi 1e-20) |
| 実行時間 | 65 s vs 58-63 s |

星形成・sink は 6 step では起きない (z~105 まで)。残る既知の差: cloudy テーブル (方針 1)、`x_MIN`、初期温度の 0.5%、
未移植の `flg_thin_star_cell` / `pop2_SN_model` / `fof_find` / `flg_voffset_ini`。星形成が起きる段階 (もっと長いラン) の比較は未実施。

### 長いラン (2026-11-14、`ksref_long` / `kstest_long`、12 ランク、nstepmax=400、出力 z=30/25/20/14)

同じ設定で星形成が起きるまで走らせた (ks は -O2 + 境界チェックのバイナリなので 1 step あたり約 2 倍遅い)。

- z=30 (output_00002): セル数 2,332,672 vs 2,332,808、密度の総和 1e-4、運動量 1.5e-4、内部エネルギー 3e-4、xHII の総和 1e-3 (最大値 5%) の相対差。
  log のエネルギーは step 98 まで 1e-3 以内 (takashi の log の桁数)。初期の丸め差と温度 0.5% が構造形成とともにゆっくり広がる。
- **最初の Pop III 形成 (pop3_model=6) は両者で同じイベント**: 同じ fine step 377 (coarse step 134、a=0.04562 / 0.04563、z≈21)、
  同じ候補セル (位置 (0.47503, 0.49184, 0.49693)、nH=11.1 cm⁻³、T=1030 K、yHI 0.999、yHII 1.2e-4、yH₂ 4.0e-4、Z=0、
  rank 8 / level 16)。雲スケールの量: ir_cloud=56、R_cloud 20294 vs 20295 pc、M_cloud 2.5689e9 Msun、
  Ṁ_cloud 7.1405e-2 vs 7.1411e-2 Msun/yr、M_pop3 2412.9 vs 2413.0 Msun (5e-5)。8 連星の質量 (1210, 605, 302, 151, 76, 38, 19, 9.5 Msun)
  も同一。takashi 側は ugvtk (level 6-16 の 11 本) と `SFdata/pop3prof.378.1.dat` も出た。
- Pop III の寿命末 (303 Msun → 2.15 Myr、1211/606 Msun → 2.4 Myr; M>260 なので直接崩壊、E_SN=0) と残骸 BH の sink 化
  (303, 1211, 606 Msun の sink → 合体で 1 個 2120 Msun、Bondi 降着 2.4e-9 Msun/yr、L=1.3e37 erg/s) が takashi 側で起きた
  (fine step 449-471)。ks 側は到達待ち (下に追記)。
- **見つけたバグ (takashi の FIRE-3 feedback、merge 前から)**: `feedback.f90` の `t_age_Gyr` と `dt_Gyr` が `scale_t` で換算されて
  いたが、`use_proper_time` のとき (cosmo + rt_star) は texp と `dteff*aexp**2` の単位が `scale_t/aexp**2` (= 1/H0) なので
  年齢と dt が a² 倍小さくなる (z=20 で 2e-3 倍) → cosmo では Pop II の SN が永久に起きない。DICE (aexp=1) では無害。
  `tsec` (= scale_t または scale_t/aexp**2) に統一して修正 (2026-11-14)。Pop III の寿命判定 (`t0_3`) は元から正しく、
  log の `t_age[Myr]` だけが a² 倍小さく表示されていた。- **FIRE-3 の年齢の評価点を ks 式に変更** (2026-11-14、ユーザー判断): ks の `fire3_yields` 経路は SN 率・質量放出率を step 中央の年齢
  (`t_age - 0.5*dteff`、midpoint 則) で評価し、step の途中で生まれた星は年齢 0 からその時点までを積分する。takashi は step 終端の
  年齢で dteff 全体に掛けていた (端点則、生まれる前の時間も数える)。`feedback.f90` を ks 式にした。非 cosmo (DICE) でも
  結果が変わる (bit 不一致、SN の立ち上がりが半 step 遅く・少なくなる) ので、DICE の回帰 (fs3) の差を下に記録。  DICE での効き方 (フィードバック段階の試験 fb5 → fb6、`dice_rt_amr` 系 14 step、sf_fg、Pop II 約 5500 個、SN_reservoir):
  92 ファイル中 72 が変わり、SN 放出 9 件 1.92e54 erg → 8 件 1.76e54 erg (-8%)、新星 5483 → 5534、output_00004 の
  内部エネルギー総和 -2.8%、密度・金属の総和 1e-4。sink 試験 (fs2/fs3: Pop II 粒子なし) は bit 一致のまま。
  a² の修正だけでは DICE は bit 一致 (fs3 で確認)。- **ks 側も同じ fine step 449 で同じ星 (idp 3、302.81 Msun) が直接崩壊し、同じ位置・速度の sink 1 (303 Msun) ができた**
  (ks の log は t_age=2.14 Myr; takashi の表示は上の a² の件で 4.6e-3 と出ていた、修正済み)。sink の質量は 302.8118 vs 302.8277
  (5e-5; 形成時の M_pop3 の差)、降着の周辺量 ninf 0.658 vs 0.666 cm⁻³、Ṁ_now 4.9e-11 vs 5.1e-11 Msun/yr (4%、局所の
  ガス状態の差)、logsink1 の 1-2 行目も同程度で一致。- **続く崩壊・sink の合体も一致**: idp 1 (1211 Msun) と idp 2 (606 Msun) の直接崩壊は ks が fine step 458/459、takashi が 459/459
  (寿命 2.43 Myr の境界が 1 step 違い)。sink 2 (1211 Msun)・sink 3 (606 Msun) の質量・位置・速度は同一。3 個は 2 回の合体で
  1 個 (ks 2119.68 Msun、takashi 2119.79 Msun) になり、合体時の距離は ks 5.8e-3 / 4.0e-3 pc、takashi 1.28e-2 / 2.6e-3 pc
  (r_merge=138.6 pc に対して十分小さい)。合体直後の Bondi 降着率は ks 2.67e-9 vs takashi 2.36e-9 Msun/yr (12%、局所のガスの差)。
  この後 idp 4 (151 Msun) の PISN (1.44e52 erg、Mej=151 Msun、Z_ej=0.47) が takashi 側で fine step 475 に起きた (ks は未到達の
  うちにランを止めた: ks バイナリは -O2 + 境界チェックで 2 倍以上遅く、PISN 後は dt が 8e-5 に落ちて進まなくなるため)。
  結論: 星形成 → 寿命末 → BH sink → 合体 → 降着の連鎖が ks と同じイベント・同じ数値 (1e-4 以下、降着率は 1e-1) で再現する。

### PISN 以後の比較 (2026-11-15、`ksref_pisn` / `kstest_pisn`、12 ランク、初期条件から z=20 まで)

ks 側を `-O3`・境界チェックなしで別ビルドし (`Ramses-FG/bin_fast/`、`make -f ../makefiles/Makefile_0X.ws F90="mpif90 -O3 ..."`)、
takashi 側は HEAD (FIRE-3 の年齢修正後) で、両方を最初から z=20 まで走らせた (ks 5.6 時間、takashi 3.4 時間; 出力 a=0.0323, 0.0385,
0.0462, 0.0465, 0.0470, 0.0476)。

- **イベントは z=20 まで全部一致**: Pop III 8 連星の寿命末 7 件 (直接崩壊 303/1211/606/76 Msun、PISN 151 Msun 1.436e52 erg、
  SN 38 Msun 2.739e52 erg Mej=22.9 Msun、SN 19 Msun 1e51 erg Mej=15.2 Msun) が同じ星・同じ質量・同じエネルギーで起き、fine step は
  449/458/459/475/503/522/1016 (ks) vs 449/459/459/474/502/522/1086 (takashi; 最後の差は hypernova 後の dt の履歴、年齢は 8.85 vs 8.87 Myr)。
  SN 後の dt も同じ (PISN 後 3.8e-4 vs 3.9e-4、38 Msun の SN 後 3.8e-5 vs 3.7e-5)。z=20 で Pop III 1、BH 7、Pop II 0、
  sink 2214.09 vs 2214.21 Msun。log の `Zej` は ks 23.5 / takashi 0.470 と 50 倍違うが (未確認: 表示が Zsun 単位か絶対値かの差と推定)、
  星の周り (箱の 2% 半径) のメタル質量は 71.20 → 79.69 → 82.48 Msun (ks) vs 71.20 → 79.70 → 82.49 Msun で 5e-5 一致。
- **HII 領域は takashi が小さい (ks の試験 namelist の群 5 が原因)**: 同じ球内の電離水素質量 (ks → takashi) は a=0.0462 で
  1.83e8 → 1.60e8 Msun、0.0465 で 4.15e8 → 3.42e8、0.0470 で 8.43e8 → 5.85e8、0.0476 で 8.24e8 → 5.53e8。熱エネルギーも 11-28% 小さい。
  星の光子放出量 (log) は同じ、LW 群の飛行中の光子数は 4 桁一致、電離 3 群の飛行中の光子数は星形成直後 (a=0.0457) は 0.1% 一致で、
  a=0.0460 で 4%、0.0462 で 7% takashi が少ない (z=25 の出力から出力を足して再開した `ksref_ts` / `kstest_ts`)。
  原因: ks の試験 namelist (`run_l5.sh` が作る) は `groupL0/L1` を 4 群しか書かず、`isFUVbin=.true.` の群 5 がコードの既定値
  (`groupL0=13.6`, `groupL1=0` = 上限なし) のままになって、電離光子の帯域全体をもう 1 群ぶん運んでいる (飛行中の光子数は群 2 の
  0.55-0.62 倍、星の死後は電離群と同じく吸収されて減る)。ks の本番 namelist (`ngridmax_fire3.nml` 等) は群 5 を 5.6-11.2 eV と明示。
  「ks との比較ラン」で群 5 を「エネルギー 0 の死んだ群」としたのは誤りだった (log の初期値 egy=0 は on the fly で更新される)。
- **確認**: takashi 側の群 5 を ks と同じ定義 (13.6 eV〜上限なし、`group_on(5)=1`) にして z=25 から再開 (`kstest_g5`) すると、a=0.0462 で
  電離水素 1.78e8 Msun (ks 1.83e8、差 13% → 2.7%)、HeIII 1.37e8 (ks 1.46e8)、熱エネルギー 1.83e54 erg (ks 1.86e54)、
  飛行中の光子数は電離 3 群 1.6%、群 5 1.7% の差になる。`attn_after_chem` (ks の試験も `&PHYSICS_PARAMS` で .true.、両者で有効) と `relax_x_by_ne` (takashi 独自) を takashi 側で切っても
  ほとんど変わらない (1.60e8 → 1.57e8 Msun、`kstest_noattn`)。
  `namelist/ks/l5_14_test.nml` は意図した物理 (群 5 = FUV 帯、星からは入れない) のままにし、ks の試験を数値まで再現する 4 行を注に書いた。
- **残った差の内訳 (2026-11-15、a=0.0460 から coarse step ごとに出力した `ksref_sn` / `kstest_sn_cur` / `kstest_sn_ks`)**:
  - HII 領域の残り 2.7% は `attn_after_chem` の解き方の違い。ks も試験では `attn_after_chem=.true.` (`&PHYSICS_PARAMS`)。どちらも
    吸収前の輻射場 (注入だけ足した N*) で chemistry を解くが、ks はそのあと**更新後**の電離状態で吸収を計算し収束判定を使わない。
    takashi は**更新前**の電離状態で吸収した値を取っておいて最後に採用し、光子密度・flux の 10% 則 (`rt_skip_convergence=.false.` が既定)
    で substep を刻む。`namelist/ks/*.nml` は `rt_skip_convergence` を書いていなかったので比較ランは判定ありで走っていた
    (DICE 本番 `xd_rt.nml` と `tests/dice_rt_*.nml` は `.true.`)。ks 方式に切り替えるスイッチ `attn_new_abund` (`&RT_PARAMS`、既定
    .false. で bit 一致) を入れ、`attn_new_abund=.true.` + `rt_skip_convergence=.true.` にすると a=0.0462 で電離水素 1.825e8 Msun
    (ks 1.832e8、-0.4%)、熱エネルギー +1.0%、飛行中の電離光子 -0.2% まで合う (現行設定は -2.7% / -1.2% / -1.1%)。
  - PISN の注入は一致 (直後の出力で T>1e5 K 2.54e4 vs 2.52e4 Msun、Tmax 1.31e7 vs 1.33e7 K、中心の平均温度 62.3 vs 61.4 kK)。その後
    takashi のほうが速く冷える (T>1e5 K のガスが 2-5 step 後に 3.98e4/2.56e4/5.9e3/1.8e3 vs 4.13e4/4.33e4/2.66e4/9.2e3 Msun、中心温度
    23.6 vs 26.3 kK)。吸収の方式を変えても同じなので冷却側の差 (ejecta は Z=0.47; cloudy テーブル = 方針 1 と推定、冷却率は未比較)。
  - Pop II は z=20 までどちらも形成されず、Pop II 形成の比較は未実施。

### 吸収の順序と収束判定 (2026-11-15、`tests/ifront_thick.nml`)

1 substep で a = k N* ddt、tau = k n_old ddt (k = sigma c) と書くと、chemistry の電離数は N* tau/(1+a)、除かれる光子は
takashi N* tau/(1+tau)、ks N* tau/(1+a+tau)。ks は常に電離数より少なく除き (電離して中性でなくなった原子の分を数え損ねる)、
takashi は N* = n_old (電離波面) で保存し前後で誤差の符号が変わる。どちらも tau >> 1 では電離数が吸収数の ~(1+tau) 倍になる
(吸収前の輻射場で chemistry を解くことの誤差) ので、substep の tau を抑える収束判定が効く。

光学的に厚いセルを進む R 型波面 (一様 nH=0.1、64^3、セルの tau~14、RT 1 step の tau~4、再結合時間より短い 0.4 Myr) での電離原子数:

| 方式 | rt_courant 0.8 | rt_courant 0.1 |
|---|---|---|
| takashi の順序、判定あり (既定) | 6.856e60 | 6.809e60 |
| ks の順序 (`attn_new_abund`)、判定あり | 6.925e60 | 6.841e60 |
| takashi の順序、判定なし (DICE 本番の設定) | 9.801e60 | 7.013e60 |
| ks の順序、判定なし (ks 相当) | 1.002e61 | 7.048e60 |
| 上流の既定 (`attn_after_chem=.false.`) | 6.359e60 | 6.386e60 |

順序の差は 1-2%。判定なしは通常の刻みで 4 割以上多く電離し、刻みを 1/8 にすると判定ありの値に寄る (= 光子の作りすぎ)。
上流の既定は刻みによらず 7% 少ない (cooling の substep 内の誤差で、どちらが真値に近いかはこのテストでは決まらない)。
通常の Stromgren 球 (nH=1e-3、平衡半径) はどの方式も 1% 以内で一致し、この差は見えない。DICE 本番の `rt_skip_convergence=.true.` が
R 型の過渡期にどれだけ効くか (セルの tau は本番のほうが大きい) は未評価。
- 道具: `tests/rt_photon_total.py` (rt ファイルの level 5 のセルから群ごとの飛行中の光子総数)。yt は sink のない出力で空の `sink_*.csv` に
  落ちる (`EmptyDataError`)、ks の出力は `hydro_file_descriptor.txt` が実態と違うので `yt.load(..., fields=[...])` で変数名を与える。
  ks は restart ファイルと namelist の `noutput` の一致を要求する (途中で `aout` を足すときは個数を保つ)。


## get_27_neighbors: 低い側の隣を 1 セル飛ばしていたバグ (2026-11-22、Ramses-FG 56e45b5 から)

- `star_formation.f90` の `get_27_neighbors` (1 セル版) と `get_27_neighbors_grid` (grid ベクトル版) は、対象が隣の grid にある
  ときの子セル番号を `mod` + 「負の補正」で作っていた。`mod(-1,2) = -1` に対して補正が +1 (+2, +4) しか足さないので、
  -1 の位置が子番号 0、つまり **-1 dx ではなく -2 dx のセル**になっていた。grid 内で ix (iy, iz) = 0 のセルでは、27 セルの箱が
  低い側の隣接セルを含まず、その 1 つ先を含んでいた (8 個の子セルのうち 3 方向とも正しいのは (1,1,1) だけ)。
- 修正: `child_pos = 1 + modulo(target_i,2) + 2*modulo(target_j,2) + 4*modulo(target_k,2)` (補正の 3 行は削除)。両方のルーチン。
- **影響するのは FG だけではない**: 本番の DICE (`sf_model=6`) の速度分散 `sigma2` は `get_27_neighbors_grid` の箱で作る
  (alpha_vir と Mach 数)。sf_fg では密度の極大の判定と `loc_vir` の速度分散・箱の質量 (Pop III / Pop II / SFC 候補 / marker)。
  この修正を含むバイナリは、これまでの本番ラン (`hu_lv15` など) と星形成の条件が変わる。量的な影響は未評価。
- 確認 (DICE 4 ランク、`sf_model=6`、最細 level): 27 セルの中心位置を中心セルからのずれ (dx 単位) で出す一時的な出力で、
  修正前は ind=1 (0,0,0) のセルで低い側が x, y, z とも -2.00、ind=8 (1,1,1) は ±1.00; 修正後はどちらも -1 / 0 / +1 が 9 個ずつ。
  手元の短いテスト (`sf_model=6` 6 step、FG の SFC 8 step) は完走し、出力は修正前と bit 一致 (どちらも箱の値で結果が変わる
  星形成が起きない設定で、効果は見えていない)。

### div v の符号 (局所 Pop II、pop2_model 2 / 3; 2026-11-22、同じ commit から)

- `star_formation_fg.f90` の div v は公開版の block を写したもので、ul を隣 2、ur を隣 1 から取っていた。`getnbor` は
  1,3,5 = 左 (-x,-y,-z)、2,4,6 = 右 (`getnborcells` の `ggg(ind,1)`: ind=1 (ix=0) の方向 1 は左隣の grid)。公開版は
  `(ur-ul)**2` しか使わないので害が無いが、こちらは `divv` の符号を使うので **-div v** になり、`if(divv>0) ok=.false.` が
  収束しているセルを落とし、発散しているセルを通していた。
- 修正: ul を隣 1 / 3 / 5、ur を隣 2 / 4 / 6 から取る。`star_formation.f90` (公開版の経路、`divv**2` だけ使う) は触っていない。
- 影響: pop2_model 2 / 3 の局所 Pop II だけ。SFC モデル (pop2_model 4 / 5)、SFC marker (div v の前の条件で決める)、Pop III は
  div v を使わない。
- 確認 (`tests/dice_rt_fg.nml`、pop2_model=2、4 ランク 6 step): 新しい星が 2695 → 280 個。回転円盤の初期条件で発散している
  セルのほうが多く、修正前はそちらで星ができていた。崩壊する塊での星形成は手元では試していない。

### 変換されていない隣セルと、27 セルの箱から取る div v (2026-11-22、Ramses-FG 453a780 から)

- ks の問題は 2 種類: (a) 他ランクの ghost grid、(b) 親セル (隣の grid が無いときに返る 1 level 粗いセル)。どちらも `uold(2:4)` が
  運動量のままで、速度として足されていた。**takashi では (a) は元から無い**: `star_formation` / `star_formation_fg` は level の
  virtual grid も速度に変換している (`sf - ghostzones`)。残っていたのは (b) と、物理境界の grid (`cpu_map = 0`)。
- 実装: `get_27_neighbors(ind_cell, neighbors, is_prim, ilevel)` と `get_27_neighbors_grid(…, neighbors, is_prim, ilevel)` が
  `is_prim` を返す。同じ grid 内は真、隣の grid のセルは `cpu_map(親セル) > 0` (自ランクと ghost; ks は自ランクだけ真)、親セルは偽。
  偽のセルは密度で割って速度にする。使う所: `star_formation.f90` の sf_model 6 / 7 の `sigma2`、`star_formation_fg.f90` の
  Pop III と Pop II / SFC の箱 (2 か所)、div v。
- div v (pop2_model 2 / 3): 6 面の隣を `getnbor` ではなく `get_27_neighbors` の箱 (13,15 / 11,17 / 5,23) から取る。隣の grid が無い面は、
  `getnbor` だと自セル (その面の速度差 0)、箱だと親セル (1 dx 離れているものとして使う)。それまでの条件を通ったセルでだけ計算する。
  `getnbor` は pop2_model 3 のポテンシャル極小の判定用に残る。
- 確認 (DICE 4 ランク、levelmin 6 / levelmax 8):
  - 通常の設定では反映前と bit 一致 (`sf_model=6` 6 step、`pop2_model=2` 6 step で星 280 個、SFC 8 step)。
  - 理由: `get_27_neighbors_grid` が返す隣のうち親セルは 1.8% (109264 / 6127056) あるが、密度の条件を通るセル (n_star=1e-5 でも
    112054 個) の箱には 1 個も入らない。質量で細分化する DICE では、level の端のセルは低密度で、星形成の候補は level の内側にある。
    **本番の DICE でこの修正が効く場面はほぼ無い** (効くのは同じ日の `modulo` の修正のほう)。
  - 経路の確認: `n_star=1e-9` にして全セルを評価させると、親セルを含む箱が 233136 個 (level 6 / 7 / 8)。速度分散は修正前の
    0.25–1.00 倍 (幾何平均 0.72)。修正前は親セルの運動量 (密度 ≪ 1 なのでほぼ 0) を速度として足し、箱の平均速度との差で
    分散を水増ししていた。
  - FG の zoom (細分化が密度の峰を追う設定) での効果は未確認。


## 速度分散の条件 (`loc_vir`) を外し、div v < 0 にする (2026-11-24、Ramses-FG 00ab333 から)

Ramses-FG 側の判断 (ユーザー) をそのまま移した。理由は ks 側の `NOTES_pop2_model_comparison_2026-11-11.md` (2026-11-24 の節):
27 セルの速度分散はほとんどが落下そのもので (形成地点 4 か所で線形の速度場を当てはめると約 6 割が等方的な収縮、ひずみの
固有値はすべて負、残差 1-9%)、自由落下する箱は T = 0 でも条件を通れない。DICE の本番 (`sf_model=6`) で見た結果と同じ
(`docs/sf_model6_sigma.md`)。

新しい条件 (`star_formation_fg.f90`):

| | 条件 |
|---|---|
| Pop III | セルの条件 (密度の閾値、Z、yH2 など) + 27 セルの密度の極大 (常に; 以前は `loc_vir > 0` のときだけ) + div v < 0 |
| Pop II、`pop2_model` 2 / 3 | 熱的な Jeans 閾値 + T < Tpop2 + div v < 0 (3 はポテンシャルの極小も)。`pop2_peak_only` で極大も |
| Pop II、`pop2_model` 4 / 5 (SFC) | 同上 (div v < 0 が新しく入った)。SFC の中心を極大に限るのは `pop2_peak_only` |
| Pop II、`pop2_model` 1 | 変更なし (div v なし) |
| SFC 印 (`pop2_SFC_marker`) | kind 1: セルの条件 + div v < 0 + 密度の極大 (以前は div v なし、極大、`loc_vir`)。kind 2: 局所モデルの星形成条件 (極大は不要) |

- `loc_vir` は変数ごと削除した (ユーザーの指示、同じ日; Ramses-FG は読むだけ残している)。**`loc_vir=...` を書いた namelist は
  起動時に止まる** (「Cannot match namelist object name」) ので行を消す。`namelist/ks/` の 3 本 (`l5_14_test.nml`、
  `ngridmax_fire3.nml`、`ngridmax_music.nml`) からは消してある。
- 同時に直したバグ: `ivar_refine > 0` のとき幾何条件が `d` を refinement のスカラーで上書きし、Pop II の div v が中心セルの
  重みにそれを使っていた (面の値がほぼ隣の速度になる)。`dcell = uold(ind_cell,1)` を使い、幾何条件は `d` を上書きしない
  (Pop III と Pop II)。`ivar_refine = 0` のランは変わらない。
- 移していないもの:
  - Ramses-FG 3f702b3 (テスト問題の星の二重計上): こちらは `add_star_particle_for_test` を移植していない。
  - Ramses-FG 4c997c6 (実 SFC の粒子は印を妨げない)、78a3194 (年齢出力の idp を標準出力に): 未移植。
- **印の kind 2 (Ramses-FG 14d37df、同じ日に移植)**: 局所モデルの星形成条件 (`ok(i)`) を満たすセルも印の候補にする。
  極大でないセルで先に星ができ、そのフィードバックで極大のセルが SFC の条件を満たさなくなる (星はあるのに印も入れ子の
  出力も無い) のを避けるため。同じ step の候補が `SFC_dist_min` 以内にあるときは kind 1 を優先し、次に密度の高いほうを残す。
  有効な印から `SFC_dist_min` 以内のセルは候補に載せない (`collect_active_SFC` をセルのループの前に呼ぶ)。
  標準出力の `** SFC marker **` の行末に `, kind = 1|2` が付く。SFC モデル (`pop2_model` 4 / 5) の経路は変えていない。

確認 (`tests/dice_rt_fg.nml` 相当、4 ランク、変更前 = ca634401):

| 試験 | 変更前 | 変更後 |
|---|---|---|
| 局所 Pop II (`pop2_model=2`、6 step) | 星 280 | 星 280、出力は bit 一致 (20 ファイル) |
| SFC 印 (上 + `pop2_SFC_marker`、`SFC_dist_min_pc=5d3`) | 印 1 個 | 同じセルに印 1 個 (kind 1)、出力は bit 一致 |
| SFC 印 (`SFC_dist_min_pc=1d3`) | — | 印 11 個 (kind 1 が 1 個、kind 2 が 10 個)、互いに 1 kpc 以上離れている。星 280 個は同じ |
| Pop III (`Zpop3=1d1`、`pop3_model=1`、6 step、`loc_vir` 未設定) | 846 個 (極大も速度分散も見ない) | 2 個 (極大 + div v < 0) |
| 同上、変更前に `loc_vir=1` | 2 個 (極大 + 速度分散)、変更後と同じ 2 セル | — |
| SFC (`pop2_model=4`、8 step) | `logSFC` 16 行、星 214 | `logSFC` 8 行、星 292 |

Pop III が 846 → 2 個に減るのは、極大の判定が `loc_vir > 0` のときだけだったのを常に課すようにしたため (試験の namelist は
`loc_vir` を設定していない)。`loc_vir=1` で流していたランと比べれば条件は緩む向き (速度分散の条件が div v < 0 に替わる) だが、
この試験では極大のセルが 2 つしかなく、どちらも両方の条件を通るので差は出ない。SFC は条件が変わるので結果も変わる。div v と極大の値そのものはこちらでは数値を突き合わせていない
(ks 側で解析的な速度場の箱を使い、符号が 1728/1728 セルで一致することを確認済みの同じコード)。`ivar_refine > 0` の経路は
こちらでも未確認。`tests/dice_rt_fg.nml` のコメントにある「Zpop3=1d1 で 846 個」を前提にした試験 (SN、sink、ugvtk) は
星の数が変わる。

## MUSIC 初期条件の本番 (hybrid、C_eps1 の物理) と、そのための移植 3 件 (2026-11-27)

新しい初期条件 `music_M1e10_in_M5e11` (4sigmaM10 の後継) と `music_M5e10_in_M1e12` (5sigmaM5e10 の後継;
`~/work7/rare_peak/music_ic_settings.md`) を patch/takashi で走らせるための一式。物理は ks の Pop II 比較ラン C_eps1
(`4sigmaM10/pop2_compare/C_eps1/C_eps1.sh` = 本番 `ngridmax_intel.sh` の sed + 局所 Pop II の 5 項目) にそろえた。

| ファイル | 内容 |
|---|---|
| `bin/Makefile.xd_ks_intel` | 既存のものをそのまま使う (oneAPI `mpiifx -qopenmp`、`OPENMP = 1`)。`bin_ks/` にコピーして組む |
| `namelist/ks/music_M1e10_in_M5e11.nml` | levelmin 8 / levelmax 23、初期条件 level 8-12 |
| `namelist/ks/music_M5e10_in_M1e12.nml` | levelmin 8 / levelmax 24、初期条件 level 8-13 |
| `jobfiles/xd/ks/music_M1e10_in_M5e11.sh` | 4 ノード、16 ランク/ノード x 7 スレッド = 64 ランク |
| `jobfiles/xd/ks/music_M5e10_in_M1e12.sh` | 16 ノード、同 = 256 ランク |

- ノード数は元のラン (flat 448 / 1792 ランク) と同じ。`ngridmax` はランクあたりなので flat の 7 倍 (2.8e6 / 3.15e6; 総量は同じ)。
  `npartmax` は両方 3.2e7 (flat の 7 倍は 2.24e7 / 1.4e7)。初期の粒子は 9.8e7 / 3.8e8 個で、zoom の外の境界粒子が
  2.0e7 / 2.6e7 個 (level 8 の 256^3 だけで 1.67e7) ある。コスト重みの分割では安い境界領域が少数のランクに集まるので、
  level 8 が 1 ランクに全部乗っても入る大きさにした。メモリは 1 ノード 16 ランクで 300-320 GB の見積もり (`--mem=461G`)。
- 光子群はエネルギーの低い順 (ユーザー 2026-11-27; やり直しのランなので ks と違ってよい): 1 FUV 5.6-11.2、2 LW 11.2-13.6、
  3 HI、4 HeI、5 HeII (ks は LW, HI, HeI, HeII, FUV)。`LWgroup=2`、`iPEH_group=1`、`spec2group=2,3,4,5`。群の番号は全部
  namelist で渡すのでコードは同じ。この並びでの実行は未確認 (読み込みと群の表の表示まで)。
- 分割と重力 (ユーザー 2026-11-27): `two_tier_level=-1` (Ls を自動選択)、`cg_auto=.true.`、`cg_warm_start=.false.`
  (warm start は高赤方偏移の zoom では逆効果)。`omp_cg=.true.` (ランクあたりの grid が多いので CG もスレッド化)。
- C_eps1 の 5 項目: `pop2_model=2`、`eps_star=1`、`pop2_peak_only=.false.`、`pop2_SFC_marker=.true.`、
  `SFC_output_at_given_age=.true.`。
- 出力: ks の `dtexpout_Myr=5` は未移植なので、開始 (z = 127、宇宙年齢 11.9 Myr) から 5 Myr ごとの `aout` を 233 個並べた
  (最後は `aend=0.167`; RAMSES と同じく輻射なしの平坦 ΛCDM で計算)。
- job script: restart は `output_restart/` の有無で決める (`nrestart=-1` / 0)。`run.lock/` で同じディレクトリの同時実行を
  防ぐので、何本かまとめて投入できる。初期条件は `/work/$USER/smbh/ICs/<名前>/` (その下に `ics_ramses/` があればそちら)。
- 古い `~/work7/rare_peak/music_*/ngridmax_intel.sh` と `ics_ramses/ngridmax_music.nml` (2026-11-14、flat) は使わない
  (削除済みの `loc_vir` が残っていて読めない)。

### 移植 3 件 (既定値では変更前と bit 一致)

C_eps1 の namelist で `!KS` のまま残っていた項目のうち、結果に効くものを入れた。

| 項目 | 場所 | 内容 |
|---|---|---|
| `jeans_refine_nH` (`&REFINE_PARAMS`、既定 0) | `hydro_flag.f90` の `jeans_length_refine` | Jeans 細分化は nH がこの値 [H/cc] を超えるセルだけ。あわせて ks と同じく、`ivar_refine > 0` のとき zoom スカラーが `var_cut_refine` 未満のセルは流体の細分化条件を全部落とす |
| `strict_Npart_refine`、`Mmin_strict_Npart_refine` (`&REFINE_PARAMS`、既定 `.false.` / -1) | `rho_fine.f90` と `to_dice/rho_fine.f90` の 2 箇所ずつ | `m_refine` の粒子数で baryon 粒子を質量 / `mass_sph` ではなく個数で数える。`Mmin` より軽い粒子は数えない |
| `pressure_fix_cond=3`、`pressure_fix_eta` (`&HYDRO_PARAMS`、既定 0.05) | `godunov_fine.f90` の `set_uold` | mode 1 に「`e_cons < pressure_fix_eta × E_total`」を足したもの。`e_cons <= 0` は条件によらず `e_prim` で置き換える。非有限・非正のエネルギーは `clean_stop` |

- `jeans_refine_nH` と zoom の制限は本番に必須: 無いと z = 127 の一様なガス (nH = 0.4 /cc、T ≈ 300 K、Jeans 長 ≈ 29 kpc/h
  comoving < 8 セル) が箱全体で level 11 まで細分化される。今までの比較ラン (`l5_14_test.nml`) は `jeans_refine=0` だったので
  この経路を通っていなかった。
- `no_sinkcloud_refine` は移植済みだったが `&SINK_PARAMS` で読む (ks は `&REFINE_PARAMS` 相当)。本番 namelist で `.true.` にした。
- `mhd/godunov_fine.f90` には `pressure_fix_cond=3` を入れていない。
- 確認: ks 設定 (gfortran、OpenMP なし) で `namelist/ks/l5_14_test.nml` (16 ランク 6 step) の出力 384 ファイルが変更前と bit 一致。

### C_eps1 と違うまま残るもの

`!KS` で残っている項目は ks 本番でもオフか既定値のもの (`fof_find`、`fout_*`、`met_switch`、`flg_thin_star_cell`、`nHmax_preSN`、
`fJeans_cool`、`SS_LVG`、`isXray`、`p2_no*`、`sgagn_seed`、`f_rho_add`、`flg_voffset_ini`) で、物理の差にはならない。差が残るのは
これまでの方針どおり: フィードバックの実装 (takashi の `feedbk` + FIRE-3 yields + `msp0`; ks の `pop2_SN_model='fixed'`、
`pop2_sESN` は無い)、cloudy の金属冷却テーブル (方針 1)、光電加熱の効率 (`iPEH_group`)、`rt_c_fraction` がスカラーなこと。
