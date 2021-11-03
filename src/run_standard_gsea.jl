"""
Run standard GSEA
"""
function run_standard_gsea(
    js::String,
    gm::String,
    tst::String,
    tsd::String,
    ou::String,
)::Any

    ke_ar = DictExtension.read(js)

    se_fe_ = read_set(gm, ke_ar)

    sa_va = TableAccess.read(tst)

    sc_fe_sa = TableAccess.read(tsd)

    return ke_ar, se_fe_, sa_va, sc_fe_sa

end

export run_standard_gsea
