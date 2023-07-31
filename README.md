The ✨ new ✨ gene-set-enrichment analysis (GSEA) 🧬.

💁 Join the [bioinformatics community](https://discord.gg/tKh7fguMrD) to get live help on GSEA (and everything bioinformatics) 🎪 ⌨️ 🔰 💡

## Use `gsea` command-line interface

#### Run single-sample GSEA

```bash
gsea data-rank
```

#### Run pre-rank GSEA

```bash
gsea user-rank
```

#### Run standard GSEA

```bash
gsea metric-rank
```

## Try with an example data

#### 1. Go to the directory with the example

```bash
cd data/sarcopenia

ls -l
```

#### 2. Make a directory for saving outputs

```bash
rm -rf output

mkdir output
```

#### 3. Run standard GSEA

```bash
gsea metric-rank setting.json target_x_sample_x_number.tsv feature_x_sample_x_number.tsv set_features.json output

ls -l output
```

#### 4. Look at the results

```bash
head -3 output/*.tsv

open output/plot/*.html
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
