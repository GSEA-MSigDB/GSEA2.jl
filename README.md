# GSEA.jl

The :sparkles: **new** :sparkles: Gene Set Enrichment Analysis :dna:

:construction: Under rapid development :construction:

**:information_desk_person: Join the [bioinformatics community](https://discord.gg/Q8XyvCfH) to get live help on GSEA (and everything bioinformatics) :circus_tent: :keyboard: :beginner: :bulb:**

## Use `gsea` command-line interface

#### Run single-sample GSEA

```bash
gsea si
```

#### Run pre-rank GSEA

```bash
gsea pr
```

#### Run standard GSEA

```bash
gsea st
```

#### Convert `.gct` and `.cls` to `.tsv`s

```bash
gsea convert-gct-and-cls bad.gct output/
```

#### Convert `.gmt` to `.json`

```bash
gsea convert-gmt bad.gmt good.json
```

#### Settings are just a [`.json` file](gsea_setting.json).

## Try with an example data

#### 1. Go to the directory with the example

```bash
cd test/data/sarcopenia

ls -l
```

#### 2. Make a directory for saving outputs

```bash
rm -rf output

mkdir output
```

#### 3. Convert `.gct` and `.cls` into `.tsv`

```bash
gsea convert-gct-and-cls gse111016_allsamplescounts_htseqcov1_sss_forgeo.sarcopenia.vs.normal_counts_collapsed_to_symbols.gct sarcopenia_binary.cls output

ls -l output
```

#### 4. Convert `.gmt` into `.json`

```bash
gsea convert-gmt c2.cp.wikipathways.v7.4.symbols.gmt output/set_to_genes.json

ls -l output

head output/set_to_genes.json
```

#### 5. Run standard GSEA

```bash
gsea run-standard-gsea ../../../gsea_setting.json output/set_to_genes.json output/target_by_sample.tsv output/gene_by_sample.tsv output

ls -l output

head output/set_by_statistic.tsv
```

## Install

1. Download a release and decompress it.

2. Add `build/gsea/bin` to the path.

3. Test installation

Start a new bash just in case.

```bash
gsea --help
```

:tada:

If there is no release matching desired machine or if installation fails, try building.

## Build

#### 1. Download this repository

```bash
git clone https://github.com/KwatMDPhD/GSEA.jl &&

cd GSEA.jl &&

julia --project --eval "using Pkg; Pkg.instantiate()"
```

#### 2. Test before building

```bash
julia --project --eval "using Pkg; Pkg.test()"
```

#### 3. Build a personal or transferable binary

###### Build a personal one

```bash
julia --project deps/build.jl
```

:point_up: installs `gsea` into `~/.julia/bin`.

If not already, add this `bin` to the path by adding :point_down: to the profile (`~/.zbashrc`, `~/.babashrc`, ...)

```bash
PATH=~/.julia/bin:$PATH
```

###### Build a transferable one

```bash
julia --project deps/build.jl app tarball
```

:point_up: makes `build` and `gsea-application-N.N.N-julia-N.N.N-MACHINE-xNN.tar.gz`.

Add `build/gsea/bin` to the path.

#### 4. Test build

Start a new bash just in case.

```bash
gsea --help
```

:tada:

---

## Howdy :wave: :cowboy_hat_face:

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

_Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
