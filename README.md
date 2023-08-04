The official command-line program for the gene-set-enrichment analysis (GSEA) 🏔️.

## Use `gsea` command-line interface

#### Run metric-rank (standard) GSEA

```bash
gsea metric-rank
```

#### Run user-rank (pre-rank) GSEA

```bash
gsea user-rank
```

#### Run data-rank (single-sample) GSEA

```bash
gsea data-rank
```

#### Convert `.cls` and `.gct` to `.tsv`s

```bash
gsea convert-cls-gct
```

#### Convert one or more `.gmt`s to a `.json`

```bash
gsea convert-gmt
```

## Try with an example data

#### 1. Go to the example directory

```bash
cd example/sarcopenia

ls -l
```

#### 2. Make a directory for saving outputs

```bash
mkdir ~/Downloads/gsea
```

#### 3. Run GSEA

```bash
gsea metric-rank \
    ~/Downloads/gsea \
    target_x_sample_x_number.tsv \
    feature_x_sample_x_number.tsv \
    set_features.json \
    --number-of-permutations 10 \
    --more-sets-to-plot "WP_DNA_MISMATCH_REPAIR WP_CELL_CYCLE ALIEN"
```

#### 4. Look at the results

```bash
cd ~/Downloads/gsea

ls -l

head -3 *.tsv

open *.html
```

## Install

1. Download the latest [release](https://github.com/KwatMDPhD/GSEA.jl/releases/latest) and extract it.

2. Add `gsea/bin` to the path.

3. Test installation

We plan to sign `gsea` in the near future. Meanwhile, enable 3rd-party apps on your macOS.

```bash
gsea --help
```

🎉

## Build

If installation is unavailable or fails, try building.

#### 1. Download this repository

```bash
git clone https://github.com/KwatMDPhD/GSEA.jl
```

#### 2. Instantiate

```bash
cd GSEA.jl &&

julia --project --eval "using Pkg; Pkg.instantiate()"
```

#### 3. Build

```bash
julia --project deps/build.jl app tarball
```

☝️ makes `build` (and `gsea-*.tar.gz`).

Add `build/gsea/bin` to the path.

#### 4. Test

```bash
gsea --help
```

🎊

---

Powered by https://github.com/KwatMDPhD/Kata.jl 🥋.
