# GSEA.jl

The :sparkles: **new** :sparkles: Gene Set Enrichment Analysis :dna:

:construction: Under rapid development :construction:

**:information_desk_person: Join the [bioinformatics community](https://discord.gg/Q8XyvCfH) to get live help on GSEA (and everything bioinformatics) :circus_tent: :keyboard: :beginner: :bulb:**

## Use `gsea` command-line interface

#### Run single-sample GSEA

```bash
gsea single-sample
```

#### Run pre-rank GSEA

```bash
gsea pre-rank
```

#### Run standard GSEA

```bash
gsea standard
```

## Try with an example data

#### 1. Go to the directory with the example

```bash
cd test/sarcopenia

ls -l
```

#### 2. Make a directory for saving outputs

```bash
rm -rf output

mkdir output
```

#### 3. Run standard GSEA

```bash
gsea standard ../../settings.json set_to_genes.json target_by_sample.tsv gene_by_sample.tsv output

ls -l output

head output/set_by_statistic.tsv
```

## Settings are just a [`.json` file](settings.json)

#### `metric`
#### `remove_gene_set_genes`
#### `minimum_gene_set_size`
#### `maximum_gene_set_size`
#### `weight`
#### `algorithm`
#### `number_of_jobs`
#### `permutation`
#### `random_seed`
#### `number_of_permutations`
#### `number_of_extreme_gene_sets_to_plot`
#### `gene_sets_to_plot`

## Install

1. Download the latest [release](https://github.com/KwatMDPhD/GSEA.jl/releases/latest) and decompress it.

2. Add `gsea/bin` to the path.

3. Test installation

```bash
gsea --help
```

:tada:

## Build

If there is no release matching desired machine or if installation fails, try building.

#### 1. Download this repository

```bash
git clone https://github.com/KwatMDPhD/GSEA.jl
```

#### 2. Download dependencies

```bash
cd GSEA.jl &&

julia --project --eval "using Pkg; Pkg.instantiate()"
```

#### 2. Test before building

```bash
julia --project --eval "using Pkg; Pkg.test()"
```

#### 3. Build

```bash
julia --project deps/build.jl app tarball
```

:point_up: makes `build` and `gsea-application-N.N.N-julia-N.N.N-MACHINE-xNN.tar.gz`.

Add `build/gsea/bin` to the path.

#### 4. Test build

```bash
gsea --help
```

:tada:

---

## Howdy :wave: :cowboy_hat_face:

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

_Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
