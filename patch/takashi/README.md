# patch/takashi

岡本研フォークの変更はすべてこのディレクトリ以下にある。作業の約束はリポジトリ直下の `CLAUDE.md`、
設計メモは `docs/`、やり残しと既知の問題は `TODO.md`。

## 使用上の注意

### `IOGROUPSIZEREP>0` (2 段の出力ディレクトリ) は非推奨

`IOGROUPSIZEREP` は `amr_parameters.f90` の中の値で、namelist では変えられない。既定の 0 (1 段、
`output_XXXXX/` と `output_restart/` にすべてのファイルを置く) のまま使うこと。

N>0 にすると出力が N ランクごとの `output_XXXXX/group_YYYYY/` に分かれるが、次の理由で勧めない。

- 上流では 2019 年からこの機能が壊れていた (ランク 1 のグループのディレクトリしか作られず、restart 時の
  info ファイルの確認も 1 段しか見ない)。このフォークでは 2026-11-22 に直したが (`ff4f35f2`)、確認したのは
  gfortran ビルドの DICE RT 4 ランクのテストだけである。
- Intel ビルド (`-DNOSYSTEM`、`PXFMKDIR` を使う側) では一度も実行していない。
- sink や ugvtk を使うランでの 2 段出力は確認していない。
- Mera.jl (`repartition` を含む) と yt が 2 段の output を読めるかは確認していない。並列度を変えた restart は
  Mera.jl に頼っているので、読めなければその手段を失う。

詳細は `TODO.md` の「4. 運用上のメモ」。

### `sf_model=6` の `sf_b_turb` の既定値が変わった (2026-11-24)

既定値を 0.4 から 0.4√3 = 0.6928 に変えた (1D の Mach 数に合わせるため)。`sf_b_turb` を書いていない namelist は、
このコミット以降のバイナリで星形成率が大きく変わる (本番の output_00038 の状態で瞬間の星形成率が約 47 倍)。以前のランを
続ける・再現するときは `&SF_PARAMS` に `sf_b_turb=0.40000000596046448` を書くこと。

同時に足した `sf_sigma_noconv` (σ から収縮流を除く) と `sf_alpha_meandens` (α_vir に 27 セルの平均密度) は既定で無効。
詳細は `docs/sf_model6_sigma.md`。
