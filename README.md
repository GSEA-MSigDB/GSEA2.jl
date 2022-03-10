# GSEA.jl

The :sparkles: **new** :sparkles: Gene Set Enrichment Analysis :dna:

:construction: Under rapid development :construction:

**:information_desk_person: Join the [bioinformatics community](https://discord.gg/tKh7fguMrD) to get live help on GSEA (and everything bioinformatics) :circus_tent: :keyboard: :beginner: :bulb:**

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
gsea standard setting.json set_genes.json target_x_sample_x_number.tsv gene_x_sample_x_score.tsv output

ls -l output

head output/set_x_statistic_x_number.tsv
```

## Settings are just a [`.json` file](setting.json)

- `metric` for ranking genes (for `standard`)

  `signal_to_noise_ratio` | `difference_of_median` | `pearson_correlation` | `information_coefficient` | `cosine_distance`

- `remove_gene_set_genes` that are not in the gene-x-sample-x-score genes

  `false` | `true`

- `minimum_gene_set_size` that removes sets smaller than this

  Integer

- `maximum_gene_set_size` that removes sets bigger than this

  Integer

- `power` to raise the gene scores

  Number

- `algorithm` for computing enrichment

  `kolmogorov_smirnov` | `jensen_shannon`

- `number_of_jobs`

  Integer

- `permutation` for computing significance

  `sample` (for `standard`) | `set` (for `standard` and `pre-rank`)

- `random_seed` from which to generate randomness for reproducibility (for `standard` and `pre-rank`)

  Integer

- `number_of_permutations` (for `standard` and `pre-rank`)

  Integer

- `number_of_extreme_gene_sets_to_plot` (for `standard` and `pre-rank`)

  Integer

- `gene_sets_to_plot` no matter what (for `standard` and `pre-rank`)

  List of string gene-set names

## Install

1. Download the latest [release](https://github.com/KwatMDPhD/GSEA.jl/releases/latest) and decompress it.

2. Add `gsea/bin` to the path.

3. Test installation

We plan to sign `gsea` in the near future. Meanwhile, enable 3rd-party apps on your macOS.

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

## :wave: :cowboy_hat_face: Howdy

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

_Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
