# GSEA.jl

The official code for Gene Set Enrichment Analysis :dna: :mountain:

![gsea2.gif](media/gsea2.gif)

![cli_example.gif](media/cli_example.gif)

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

## Use in julia

```jl
using GSEA
```

### Run single-sample GSEA

```jl
run_single_sample_gsea()
```

### Run pre-rank GSEA

```jl
run_pre_pank_gsea()
```

### Run standard GSEA

```jl
run_standard_gsea()
```

## Install

```sh
git clone https://github.com/KwatMDPhD/GSEA.jl &&

cd GSEA.jl &&

julia --project deps/build.jl
```

:point_up: commands install `gsea` into `~/.julia/bin`.

If not already, add this `bin` to the path by adding :point_down: to the profile (`~/.zshrc`, `~/.bashrc`, ...):

```sh
PATH=~/.julia/bin:$PATH
```

Start a new shell to load the updated profile.

Test installation:

```sh
gsea -h
```

:tada:

---

## Howdy :wave: :cowboy_hat_face:

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/KwatMDPhD/GSEA.jl/issues/new/choose).

---

_Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
