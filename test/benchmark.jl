using Test: @test

using BioLab

using GSEA

# ---- #

const DI = joinpath(dirname(@__DIR__), "benchmark", "outer_loop")

const IN = BioLab.Path.remake_directory(joinpath(DI, "input"))

# ---- #

const TST, TSF = GSEA.convert_cls_gct(
    joinpath(IN, "target_x_sample_x_number.tsv"),
    joinpath(IN, "feature_x_sample_x_number.tsv"),
    joinpath(DI, "Coller_et_al_phen.cls"),
    joinpath(DI, "Coller_et_al_gene_exp_preproc.gct"),
)

const JS = GSEA.convert_gmt(
    joinpath(IN, "set_features.json"),
    joinpath(DI, "h.all.v2022.1.Hs.symbols.gmt"),
    joinpath(DI, "c2.all.v2022.1.Hs.symbols.gmt"),
)

# ---- #

const OU = GSEA.metric_rank(
    BioLab.Path.remake_directory(joinpath(DI, "output")),
    TST,
    TSF,
    JS;
    normalization_dimension = 1,
    normalization_standard_deviation = 3.0,
    permutation = joinpath(DI, "KS_SUP example mean scaling_rand_perm_gene_scores.txt"),
)

# ---- #

const PME = sort!(
    BioLab.DataFrame.read(
        joinpath(DI, "KS_SUP example mean scaling_gene_selection_scores.txt");
        select = [1, 2],
    ),
)

const JME = sort!(BioLab.DataFrame.read(joinpath(OU, "feature_x_metric_x_score.tsv")))

@test size(PME, 1) == size(JME, 1)

for id in 1:size(PME, 1)

    @test PME[id, 1] == JME[id, 1]

    @test isapprox(PME[id, 2], JME[id, 2]; atol = 1e-6)

end

# ---- #

const SE_ = [1, 5, 4, 6, 7]

# ---- #

const PEN = sort!(
    BioLab.DataFrame.read(
        joinpath(DI, "KS_SUP example mean scaling_GSEA_results_table.txt");
        select = SE_,
    ),
)

const JEN = sort!(BioLab.DataFrame.read(joinpath(OU, "set_x_statistic_x_number.tsv")))

for id in 1:size(PEN, 1)

    @test PEN[id, 1] == JEN[id, 1]

    @test isapprox(PEN[id, 3], JEN[id, 2]; atol = 1e-3)

    @test isapprox(PEN[id, 2], JEN[id, 3]; atol = 1e-2)

    @test isapprox(PEN[id, 4], JEN[id, 4]; atol = 1e-2)

    @test isapprox(PEN[id, 5], JEN[id, 5]; atol = 1e-2)

end

# ---- #

const DIB = joinpath(dirname(@__DIR__), "benchmark")

const DIJ = joinpath(DIB, "json_files")

const DIT = joinpath(DIB, "Datasets_and_Phenotypes")

const DIS = joinpath(DIB, "Gene_Sets_Collections")

const DIR = joinpath(DIB, "results_sets")

# ---- #

function parse_float(fl)

    if fl isa AbstractString

        fl = parse(Float64, lstrip(fl, '≤'))

    end

    fl

end

# ---- #

const AL_ = ("ks", "kli", "kli", "kliom", "kliop")

for (id, js) in enumerate(BioLab.Path.read(DIJ))

    ke_va = BioLab.Dict.read(joinpath(DIJ, js))[chop(js; tail = 5)]

    # TODO: Debug.
    if js in ("GSE121051 LPS vs. CNTRL.json", "NRF2_mouse_model.json")

        continue

    end

    @info "$id $js"

    di = BioLab.Path.remake_directory(joinpath(DIB, BioLab.Path.clean(chop(js; tail = 3))))

    tst = joinpath(di, "target_x_sample_x_number.tsv")

    tsf = joinpath(di, "feature_x_sample_x_number.tsv")

    js = joinpath(di, "set_features.json")

    @info "Converting"

    GSEA.convert_cls_gct(tst, tsf, joinpath(DIT, ke_va["cls"]), joinpath(DIT, ke_va["ds"]))

    GSEA.convert_gmt(js, (joinpath(DIS, gm) for gm in ke_va["gene_sets_collections"])...)

    dip = joinpath(DIR, basename(ke_va["results_directory"]))

    for (al, pr) in zip(AL_, ke_va["results_files_prefix"])

        @info "Comparing $al"

        dio = BioLab.Path.remake_directory(joinpath(di, "output_$al"))

        GSEA.metric_rank(dio, tst, tsf, js; algorithm = al)

        @info "Comparing metrics"

        tx = joinpath(dip, "$(pr)_gene_selection_scores.txt")

        pda = sort!(BioLab.DataFrame.read(tx; select = [1, 2]))

        jda = sort!(BioLab.DataFrame.read(joinpath(dio, "feature_x_metric_x_score.tsv")))

        @test size(pda, 1) == size(jda, 1)

        for id in 1:size(pda, 1)

            @test pda[id, 1] == jda[id, 1]

            # TODO: Check directionality.
            @test isapprox(abs(pda[id, 2]), abs(jda[id, 2]); atol = 1e-5)

        end

        @info "Comparing enrichments"

        GSEA.user_rank(
            dio,
            tx,
            js;
            algorithm = al,
            permutation = joinpath(dip, "$(pr)_rand_perm_gene_scores.txt"),
        )
        # TODO
        continue

        pda = sort!(
            BioLab.DataFrame.read(joinpath(dip, "$(pr)_GSEA_results_table.txt"); select = SE_),
        )

        jda = sort!(BioLab.DataFrame.read(joinpath(dio, "set_x_statistic_x_number.tsv")))

        @test size(pda, 1) == size(jda, 1)

        for id in 1:size(pda, 1)

            @test pda[id, 1] == jda[id, 1]

            # TODO: Check directionality.
            @test isapprox(abs(pda[id, 3]), abs(jda[id, 2]); atol = 1e-3)

            @test isapprox(pda[id, 2], jda[id, 3]; atol = 1e-2)

            @test isapprox(parse_float(pda[id, 4]), jda[id, 4]; atol = 1e-2)

            @test isapprox(parse_float(pda[id, 5]), jda[id, 5]; atol = 1e-1)

        end

    end

end
