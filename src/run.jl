using Kwat.feature_set_enrichment: score_set
using Kwat.gmt: read as gmt_read
using Kwat.json: read as json_read
using Kwat.table: read as table_read, write


"""
Run single-sample GSEA.
"""
@cast function run(
    ke_ar::String,
    sc_fe_sa::String,
    se_fe_::String,
    ou::String,
)::ke_ar = json_read(ke_ar)

    sc_fe_sa = table_read(sc_fe_sa)

    se_fe_ = select(gmt_read(se_fe_), ke_ar["mi"], ke_ar["ma"])

    write(ou, score_set(sc_fe_sa, se_fe_; ke_ar...))

    return ke_ar, sc_fe_sa, se_fe_

end

"""
Run prerank GSEA.
"""
@cast function run(
    ke_ar::String,
    fe_va::String,
    se_fe_::String,
    ou::String,
)::ke_ar = json_read(ke_ar)

    fe_va = table_read(fe_va)

    se_fe_ = select(gmt_read(se_fe_), ke_ar["mi"], ke_ar["ma"])

    return ke_ar, fe_va, se_fe_

end

"""
Run GSEA.
"""
@cast function run(
    ke_ar::String,
    sa_va::String,
    se_fe_::String,
    ou::String,
)::ke_ar = json_read(ke_ar)

    sa_va = table_read(sa_va)

    se_fe_ = select(gmt_read(se_fe_), ke_ar["mi"], ke_ar["ma"])

    return ke_ar, sa_va, se_fe_

end
