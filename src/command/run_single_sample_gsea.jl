"""
Run single-sample GSEA

# Arguments

  - `js`:
  - `se`:
  - `ts`:
  - `ou`: output directory
"""
@cast function run_single_sample_gsea(js, se, ts, ou)

    ke_ar = dict_read(js)

    se_fe_ = read_set(se, ke_ar)

    sc_fe_sa = table_read(ts)

    en_se_sa = score_set(sc_fe_sa, se_fe_; symbolize_key(ke_ar)...)

    table_write(ou, en_se_sa)

end
