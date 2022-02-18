# GSEA.jl

The :sparkles: **new** :sparkles: Gene Set Enrichment Analysis :dna:

:construction: Under rapid development :construction:

**:information_desk_person: Join the [bioinformatics community](https://discord.gg/E49y3yhG) to get live help on GSEA (and everything bioinformatics) :circus_tent: :keyboard: :beginner: :bulb:**

## Use `gsea` command-line interface

#### Run single-sample GSEA

```
gsea single-sample
```

#### Run pre-rank GSEA

```
gsea pre-rank
```

#### Run standard GSEA

```
gsea standard
```

## Try with an example data

#### 1. Go to the directory with the example

```
cd test/sarcopenia

ls -l
```

#### 2. Make a directory for saving outputs

```
rm -rf output

mkdir output
```

#### 3. Run standard GSEA

```
gsea standard settings.json set_genes.json number.target_x_sample.tsv score.gene_x_sample.tsv output

ls -l output

head output/float.set_x_statistic.tsv
```

## Settings are just a [`.json` file](settings.json)

(`settings.json` descriptions are coming soon...)

## Install

1. Download the latest [release](https://github.com/KwatMDPhD/GSEA.jl/releases/latest) and decompress it.

2. Add `gsea/bin` to the path.

3. Test installation

```
gsea --help
```

:tada:

## Build

If there is no release matching desired machine or if installation fails, try building.

#### 1. Download this repository

```
git clone https://github.com/KwatMDPhD/GSEA.jl
```

#### 2. Download dependencies

```
cd GSEA.jl &&

julia --project --eval "using Pkg; Pkg.instantiate()"
```

#### 2. Test before building

```
julia --project --eval "using Pkg; Pkg.test()"
```

#### 3. Build

```
julia --project deps/build.jl app tarball
```

:point_up: makes `build` and `gsea-application-N.N.N-julia-N.N.N-MACHINE-xNN.tar.gz`.

Add `build/gsea/bin` to the path.

#### 4. Test build

```
gsea --help
```

:tada:

---

## :wave: :cowboy_hat_face: Howdy

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

_Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
