using Kwat.dictionary: make_keyword_argument
using Kwat.feature_set_enrichment: score_set
using Kwat.gmt: read as gmt_read
using Kwat.json: read as json_read
using Kwat.table: read as table_read, write

function _read_set(
    se_fe_::String,
    ke_ar::Dict{String, Any},
)::Dict{String, Vector{String}}

    return select(gmt_read(se_fe_), pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

end


"""
Run single-sample GSEA.
"""
@cast function run_si(
    ke_ar::String,
    se_fe_::String,
    sc_fe_sa::String,
    ou::String,
)::Nothing

    ke_ar = json_read(ke_ar)

    se_fe_ = _read_set(se_fe_, ke_ar)

    sc_fe_sa = table_read(sc_fe_sa)

    en_se_sa = score_set(sc_fe_sa, se_fe_; make_keyword_argument(ke_ar)...)

    write(ou, en_se_sa)

    return

end

"""
Run prerank GSEA.
"""
@cast function run_pr(
    ke_ar::String,
    se_fe_::String,
    fe_va::String,
    ou::String,
)::Any

    ke_ar = json_read(ke_ar)

    se_fe_ = _read_set(se_fe_, ke_ar)

    fe_sc = table_read(fe_va)

    return ke_ar, se_fe_, fe_sc

end

"""
Run GSEA.
"""
@cast function run_st(
    ke_ar::String,
    se_fe_::String,
    sa_va::String,
    ou::String,
)::Any

    ke_ar = json_read(ke_ar)

    se_fe_ = _read_set(se_fe_, ke_ar)

    sa_va = table_read(sa_va)

    return ke_ar, se_fe_, sa_va

end
