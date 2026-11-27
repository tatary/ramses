# 星風エネルギーをまとめて球に入れる (`E_wind_batch`、2026-11-26)

## 動機

DICE のやり直しラン (`hu_lv15_lsig/hu_lv15`、dx = 3 pc) で、fine step の dt (中央値 0.3 kyr) を決めていたのは
SN でも RT でもなく、若い星の星風で温まったガスだった。

- `STELLAR_POPULATION_MASS=1` では 15 Msun の星粒子の 15% が `msp0 = 100 Msun` になり、その粒子だけが 44 Myr
  (`agemax`) まで星風のエネルギーを `msp0/mp0` (= 6.7) 倍して出す。質量放出は倍しない (意図した設計: `msp0 > 0` の
  粒子にだけ質量放出を負わせると星の質量が負になり得る)。
- エネルギーは毎 step、星のいる 1 セルに熱として入る。星の周り 17 pc 以内に nH ~ 1e-4、T/μ = 1e8-5e8 K のガス
  (合計 0.1 Msun 未満) ができ、その Courant 条件が levelmax の dt を決める。星形成の開始 (90 Myr) と同時に dt が
  5 kyr から 0.3 kyr に落ち、最初の SN (101 Myr) を待たない。
- 星風の周り 20 pc を除いて計算した dt は 100-130 Myr で 5 kyr、210 Myr で 1.3 kyr (RT の上限は 2.65 kyr)。
- `msp0 = 100 Msun` の粒子 1 個の星風エネルギーは 44 Myr で 2.6e49 erg (Z = 0.1 Zsun; 99% が最初の 10 Myr)、
  同じ星団の SN 約 1 発 (1e51 erg) の 2.5%。

## 方式

`&FEEDBACK_PARAMS` の `E_wind_batch` [erg] (既定 0 = 無効、今までどおり) を正にすると、

1. `msp0/mp0` 倍される星風エネルギー (年齢 `agemax` 以下) は星のセルに入れない。
2. 粒子が生まれてから出した累積エネルギー `msp0 × W(年齢, Z)` が `E_wind_batch` の倍数を越えた step で、越えた
   回数 × `E_wind_batch` を 1 つのパケットにする。`agemax` に達した step で端数を出す。
3. パケットは粒子の位置を中心とする半径 `Ncell_wind` × (星のいるセルの幅) の球 (セル中心が球に入る葉セル) に
   体積重みで入れる (星風は等方的なので体積重み; ユーザー 2026-11-26)。既定 `Ncell_wind = 4`。経路は SN reservoir
   (`sn_pending`) で、その level が次に動き出すとき (`prepare_sn_reservoir`) に `uold` に足される。半径の探索はしない。

質量・金属の放出、`agemax` より古い星の星風 (倍率なし、全粒子)、`msp0 = 0` の粒子は変えていない。

### 状態を持たない

`W(年齢, Z) = ∫ ½ mdot v² f da` [(km/s)²] は `get_wind_mass_loss_rate` の式の積分で、`f = mp/mp0` は残っている質量の
割合 (`d ln f/da = −(mdot + CC の率 × 放出質量)`; feedbk の質量放出は SN も期待値で連続的に引く。`agemax` までは
SN Ia はない)。年齢と金属量だけの関数なので、`module wind_batch` (`feedback.f90`) が初回に表 (log Z 65 点 ×
年齢 4401 点、双線形補間) を作る。パケットを出す step は `tp`, `zp`, `msp0` だけで決まり、粒子の属性も restart
ファイルの追加もない。放出待ちのエネルギーは `sn_pending` として既に restart に入っている。

- 今までの方式との差: 実際の粒子質量は dt の刻みで表の `f` とわずかにずれるので、エネルギーの総量は刻み誤差の
  分だけ違う。
- step の境界: step の始めの年齢を `a1 − dteff` で作るので、前の step の終わりの年齢と丸めの分だけずれる。累積
  エネルギーがちょうど倍数に乗る (確率 ~1e-16) とパケットが 1 個増減しうる。
- `agemax` をまたぐ step: ホストセルに入れる側は mid-step の年齢で切り替わる (元のコードのまま) ので、半 step 分の
  倍率なしの星風が重複するか抜ける (44 Myr の星風エネルギーは無視できる)。

## 条件

- `SN_reservoir=.true.` が必要 (`sn_pending` の経路を使う)。`STELLAR_POPULATION_MASS=1` のビルドが必要。
- ログ: `WIND_BATCH coarse/fine/level=... E_erg/r_pc/M_sphere_Msun=...` (パケットごと、rank 1)。
  `SN reservoir: released` の行は SN と星風の合計になる。

## 効かなくなること

球 (4 dx = 12 pc、nH 1-100 なら 200-2e4 Msun) に 1e48 erg を入れても数百 K 程度しか上がらず、すぐ冷える。星風は
実質的に効かなくなり、SN の前に nH ~ 1e-4 の空洞を掘る働きもなくなる。SN の空洞 (nH ~ 1e-4) の中にいる若い星の
パケットは今までどおり高温になりうる。

## 確認 (2026-11-26)

`run/fb6` の設定 (DICE dw1、levelmax 8、`sf_fg`、`SN_reservoir`、4 ランク 14 step) で

- `E_wind_batch = 0`: 変更前 (03ffcf10) と `amr/hydro/part/grav/rt` の 92 ファイルが bit 一致。
- `E_wind_batch = 1d48`: パケット 6243 個、合計 4.984200e52 erg。最終出力の星 280 個について
  `Σ aint(msp0 W / 1e48) × 1e48` を表から計算した値 4.984200e52 erg と一致 (端数込みの合計は 5.004e52 erg)。
  reservoir の放出は SN 4.370e53 + 星風 4.984e52 = 4.868e53 erg。

本番の解像度 (dx = 3 pc) での dt の伸びは未確認。ランク数を変えたときの一致、restart をまたぐ一致も未確認。
