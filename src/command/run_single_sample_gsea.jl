"""
Run single-sample GSEA

# Arguments

  - `js`:
  - `se`:
  - `ts`:
  - `ou`: output directory

"""
@cast function run_single_sample_gsea(
    js::String,
    se::String,
    ts::String,
    ou::String,
)::Nothing

    ke_ar = DictExtension.read(js)

    se_fe_ = read_set(se, ke_ar)

    sc_fe_sa = TableAccess.read(ts)

    en_se_sa = FeatureSetEnrichment.score_set(
        sc_fe_sa,
        se_fe_;
        DictExtension.convert_to_keyword_argument(ke_ar)...,
    )

    TableAccess.write(ou, en_se_sa)

    return nothing

end
