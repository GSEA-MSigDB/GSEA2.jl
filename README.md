# Gene set enrichment analysis 🏔️

The Julia implementation of the next-generation GSEA is hundreds of times faster and incorporates new information-theoretic algorithms.

## Julia

```julia
]add https://github.com/GSEA-MSigDB/GSEA.jl

using GSEA

const SA = pkgdir(GSEA, "example", "sarcopenia")

GSEA.metric_rank(
    mkpath(joinpath(homedir(), "Downloads", "output")),
    joinpath(SA, "target_x_sample_x_number.tsv"),
    joinpath(SA, "feature_x_sample_x_number.tsv"),
    joinpath(SA, "set_features.json"),
    number_of_permutations = 10,
    more_sets_to_plot = ["WP_DNA_MISMATCH_REPAIR", "WP_CELL_CYCLE"],
)
```

## Command Line Interface

```bash
mkdir ~/Downloads/output

cd example/sarcopenia
gsea metric-rank \
    ~/Downloads/output \
    target_x_sample_x_number.tsv \
    feature_x_sample_x_number.tsv \
    set_features.json \
    --number-of-permutations 10 \
    --more-sets-to-plot "WP_DNA_MISMATCH_REPAIR WP_CELL_CYCLE"

cd ~/Downloads/output

head -4 *.tsv

open *.html
```

## Install

We plan to sign this app for macOS soon.
In the meantime, please [enable third-party apps](https://support.apple.com/en-us/102445#openanyway).

Download and extract the latest [release](https://github.com/GSEA-MSigDB/GSEA.jl/releases/latest).

```bash
PATH=$(pwd)/gsea/bin:$PATH

gsea --help
```

## Build

```bash
julia --project deps/build.jl app tarball

PATH=$(pwd)/build/gsea/bin:$PATH

gsea --help
```

## Contact Us

If you have any questions, issues, or concerns, please feel free to [open a GitHub issue](https://github.com/GSEA-MSigDB/GSEA2.jl/issues/new/choose).

---

Made by [Kata](https://github.com/KwatMDPhD/Kata.jl) 🥋
