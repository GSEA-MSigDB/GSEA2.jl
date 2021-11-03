# GSEA.jl

Official Gene Set Enrichment Analysis program :mountain:

## Use `gsea` command line interface

### Run standard GSEA

```sh
gsea st
```

### Run pre-rank GSEA

```sh
gsea pr
```

### Run single-sample GSEA

```sh
gsea si
```

## Use in julia

```jl
using GSEA
```

### Run standard GSEA

```jl
run_standard_gsea()
```

### Run pre-rank GSEA

```jl
run_prepank_gsea()
```

### Run single-sample GSEA

```jl
run_single_sample_gsea()
```

## Install

```sh
git clone https://github.com/KwatMDPhD/GSEA.jl &&

cd GSEA.jl &&

julia --project deps/build.jl
```

:point_up: commands install `gsea` into `~/.julia/bin`.

If not already, add this `bin` to the path by adding :point_down: to the profile (~/.zshrc, ~/.bashrc, ...):

```sh
PATH=~/.julia/bin:$PATH
```

Start a new shell to load the updated profile.

Test installation:

```sh
gsea -h
```

:tada:

## Shoutout

[Roger-luo](https://github.com/Roger-luo) made [Comonicon.jl](https://github.com/comonicon/Comonicon.jl), which powers the installation process of this project :raised_hands:

## Howdy :wave: :cowboy_hat_face:

To report a bug, request a feature, or leave a comment (about anything related to this project), just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

Made by https://github.com/KwatMDPhD/PkgRepository.jl
