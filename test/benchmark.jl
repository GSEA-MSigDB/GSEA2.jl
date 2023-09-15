using Test: @test

using BioLab

using GSEA

# ---- #

const DIB = joinpath(dirname(@__DIR__), "benchmark")

const DII = joinpath(homedir(), "Desktop", "benchmark")

const DIJ = joinpath(DII, "json_files")

const DIT = joinpath(DII, "Datasets_and_Phenotypes")

const DIS = joinpath(DII, "Gene_Sets_Collections")

const DIR = joinpath(DII, "results_sets")

# ---- #

function parse_float(fl)

    if fl isa AbstractString

        fl = parse(Float64, lstrip(fl, '≤'))

    end

    fl

end

# ---- #

function path(di)

    if RE || !isdir(di)

        BioLab.Path.remake_directory(di)

    end

    di

end

# ---- #

const RE = false

const AL_ = ["ks", "kli", "kliom", "kliop"]

# ---- #

for (id, js) in enumerate(BioLab.Path.read(DIJ))

    if js in ("GSE121051 LPS vs. CNTRL.json", "NRF2_mouse_model.json")

        continue

    end

    @info "$id $js"

    ke_va = BioLab.Dict.read(joinpath(DIJ, js))[chop(js; tail = 5)]

    dib = path(joinpath(DIB, BioLab.Path.clean(chop(js; tail = 5))))

    dii = path(joinpath(dib, "input"))

    tst = joinpath(dii, "target_x_sample_x_number.tsv")

    tsf = joinpath(dii, "feature_x_sample_x_number.tsv")

    jss = joinpath(dii, "set_features.json")

    if RE || !isfile(tst) || !isfile(tsf)

        @info "Converting `cls` and `gct`"

        GSEA.convert_cls_gct(tst, tsf, joinpath(DIT, ke_va["cls"]), joinpath(DIT, ke_va["ds"]))

    end

    if RE || !isfile(jss)

        @info "Converting `gmt`"

        GSEA.convert_gmt(jss, (joinpath(DIS, gm) for gm in ke_va["gene_sets_collections"])...)

    end

    dip = joinpath(DIR, basename(ke_va["results_directory"]))

    for (al, pr) in zip(AL_, ke_va["results_files_prefix"])

        @info "Comparing $al"

        dio = path(joinpath(dib, "output_$al"))

        feature_x_metric_x_score_tsv = joinpath(dio, "feature_x_metric_x_score.tsv")

        if RE || !isfile(feature_x_metric_x_score_tsv)

            @info "Computing metrics"

            GSEA.metric_rank(dio, tst, tsf, jss; algorithm = al, number_of_permutations = 0)

        end

        feature_x_metric_x_score = joinpath(dip, "$(pr)_gene_selection_scores.txt")

        pda = sort!(BioLab.DataFrame.read(feature_x_metric_x_score; select = [1, 2]))

        jda = sort!(BioLab.DataFrame.read(feature_x_metric_x_score_tsv))

        @test size(pda, 1) === size(jda, 1)

        for id in 1:size(pda, 1)

            @test pda[id, 1] === jda[id, 1]

            # TODO: Check directionality.
            @test isapprox(abs(pda[id, 2]), abs(jda[id, 2]); atol = 1e-5)

        end

        set_x_statistic_x_number_tsv = joinpath(dio, "set_x_statistic_x_number.tsv")

        if RE || !isfile(set_x_statistic_x_number_tsv)

            @info "Computing enrichments"

            GSEA.user_rank(
                dio,
                feature_x_metric_x_score,
                jss;
                algorithm = al,
                permutation = joinpath(dip, "$(pr)_rand_perm_gene_scores.txt"),
            )

        end

        pda = sort!(
            BioLab.DataFrame.read(
                joinpath(dip, "$(pr)_GSEA_results_table.txt");
                select = [1, 5, 4, 6, 7],
            ),
        )

        jda = sort!(BioLab.DataFrame.read(set_x_statistic_x_number_tsv))

        @test size(pda, 1) === size(jda, 1)

        for id in 1:size(pda, 1)

            @test pda[id, 1] === jda[id, 1]

            # TODO: Check directionality.
            @test isapprox(abs(pda[id, 3]), abs(jda[id, 2]); atol = 1e-3)

            @test isapprox(pda[id, 2], jda[id, 3]; atol = 1e-2)

            @test isapprox(parse_float(pda[id, 4]), jda[id, 4]; atol = 1e-2)

            @test isapprox(parse_float(pda[id, 5]), jda[id, 5]; atol = 1e-1)

        end

    end

end
