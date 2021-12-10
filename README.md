# GSEA.jl

The official Gene Set Enrichment Analysis :dna:

## Use `gsea` command line interface

### Run single-sample GSEA

```sh
gsea si
```

### Run pre-rank GSEA

```sh
gsea pr
```

### Run standard GSEA

```sh
gsea st
```

### Convert `.gct` and `.cls` to `.tsv`s

```sh
gsea convert-gct-and-cls bad.gct output/
```

### Convert `.gmt` to `.json`

```sh
gsea convert-gmt bad.gmt good.json
```

### Test sarcopenia

```sh

cd test/data/sarcopenia

gsea convert-gct-and-cls gse111016_allsamplescounts_htseqcov1_sss_forgeo.sarcopenia.vs.normal_counts_collapsed_to_symbols.gct sarcopenia_bianry.cls .

ls score.*

gsea convert-gmt c2.cp.wikipathways.v7.4.symbols.gmt set_to_genes.json

head set_to_genes.json

gsea run-standard-gsea ../setting_for_standard_gsea.json set_to_genes.json score.target_by_sample.tsv score.gene_by_sample.tsv .

head set_by_statistic.tsv
```

## Install

```sh
git clone https://github.com/KwatMDPhD/GSEA.jl &&

cd GSEA.jl &&

julia --project --eval "using Pkg; Pkg.instantiate()" &&

julia --project deps/build.jl
```

:point_up: commands install `gsea` into `~/.julia/bin`.

If not already, add this `bin` to the path by adding :point_down: to the profile (`~/.zshrc`, `~/.bashrc`, ...)

```sh
PATH=~/.julia/bin:$PATH
```

Start a new shell to load the updated profile.

Test installation

```sh
gsea -h
```

:tada:

---

## Howdy :wave: :cowboy_hat_face:

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

_Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
