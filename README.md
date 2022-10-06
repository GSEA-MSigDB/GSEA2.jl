**:information_desk_person: Join the [bioinformatics community](https://discord.gg/tKh7fguMrD) to get live help on GSEA (and everything bioinformatics) :circus_tent: :keyboard: :beginner: :bulb:**

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
cd example.sarcopenia

ls -l
```

#### 2. Make a directory for saving outputs

```bash
rm -rf output

mkdir output
```

#### 3. Run standard GSEA

```bash
gsea metric-rank setting.json target_x_sample_x_number.tsv gene_x_sample_x_score.tsv set_genes.json output

ls -l output

head -2 output/*.tsv

open output/plot/*.html
```

#### Alternatively, (instead of in command line) run this example in `julia`

```jl
using GSEA

cd("example.sarcopenia")

GSEA.metric_rank(
    "setting.json",
    "set_genes.json",
    "target_x_sample_x_number.tsv",
    "gene_x_sample_x_score.tsv",
    "output",
)
```

## Settings are just a [`.json` file](setting.json)

- `metric` for ranking genes (for `metric-rank`)

  `signal_to_noise_ratio` | `mean_difference` | `median_difference` | `pearson_correlation` | `cosine_distance` | `information_coefficient`

- `remove_gene_set_genes` that are not in the gene-x-sample-x-score genes

  `false` | `true`

- `minimum_gene_set_size` that removes sets smaller than this

  Integer

- `maximum_gene_set_size` that removes sets bigger than this

  Integer

- `exponent` to raise the gene scores

  Number

- `algorithm` for computing enrichment

  `cidac` (_cumulative information divergence with antisymmetricity and complementation_) | `ks` (_Kolmogorov Smirnov_) | `ksa` (`ks` area)

- `number_of_jobs`

  Integer

- `permutation` for computing significance

  `sample` (for `metric-rank`) | `set` (for `metric-rank` and `user-rank`)

- `random_seed` from which to generate randomness for computing significance (for `metric-rank` and `user-rank`)

  Integer

- `number_of_permutations` (for `metric-rank` and `user-rank`)

  Integer

- `number_of_extreme_gene_sets_to_plot` (for `metric-rank` and `user-rank`)

  Integer

- `gene_sets_to_plot` in addition to the extreme ones (for `metric-rank` and `user-rank`)

  List of string, gene-set names

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

If there is no release matching desired machine or installation fails, try building.

#### 1. Download this repository

```bash
git clone https://github.com/KwatMDPhD/GSEA.jl
```

#### 2. Download dependencies

```bash
cd GSEA.jl &&

julia --project --eval "using Pkg; Pkg.instantiate()"
```

#### 3. Test before building

```bash
julia --project --eval "using Pkg; Pkg.test()"
```

#### 4. Build

```bash
julia --project deps/build.jl app tarball
```

:point_up: makes `build` and `gsea-application-N.N.N-julia-N.N.N-MACHINE-xNN.tar.gz`.

Add `build/gsea/bin` to the path.

#### 5. Test build

```bash
gsea --help
```

:tada:

---

## :wave: :cowboy_hat_face: Howdy

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

Powered by https://github.com/KwatMDPhD/Kata.jl
