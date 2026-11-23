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
