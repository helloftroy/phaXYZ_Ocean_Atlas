"""Small shared statistics helpers, factored out of
phac_habitat_phylum_controlled_test.py when a second script
(plot_phac_mcl_precursor_strategy.py) needed the same stratified test.
"""
from scipy.stats import chi2


def mantel_haenszel(strata, min_stratum_n=1):
    """Cochran-Mantel-Haenszel test: combines a binary association (e.g.
    trait present/absent x in-group/out-group) across strata (e.g. one
    stratum per phylum) into one pooled odds ratio and significance test,
    controlling for whatever the stratification variable is.

    strata: iterable of (a, b, c, d) 2x2 tables, one per stratum, where
      a = trait+ in-group, b = trait- in-group,
      c = trait+ out-group, d = trait- out-group.
    Returns (OR_MH, cmh_chi2, p_value, n_strata_used).

    Formula (Mantel & Haenszel 1959); implemented directly rather than via
    statsmodels (not installed everywhere this project's scripts run) --
    short enough to audit by eye:
      OR_MH = sum(a_i*d_i/n_i) / sum(b_i*c_i/n_i)
      CMH   = (sum(a_i) - sum(E_i))^2 / sum(Var_i),  ~ chi-square, df=1
      E_i   = (a_i+b_i)(a_i+c_i)/n_i
      Var_i = (a_i+b_i)(c_i+d_i)(a_i+c_i)(b_i+d_i) / (n_i^2 (n_i-1))
    Strata where the stratification category is entirely absent from one
    side (a+b=0 or c+d=0) contribute zero to both the statistic and its
    variance, so they are skipped rather than distorting anything by
    inclusion -- mathematically equivalent to including them.
    """
    num_or = den_or = 0.0
    sum_a = sum_e = sum_var = 0.0
    n_used = 0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n <= min_stratum_n:
            continue
        if (a + b) == 0 or (c + d) == 0:
            continue
        n_used += 1
        num_or += a * d / n
        den_or += b * c / n
        e = (a + b) * (a + c) / n
        if n > 1:
            var = (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        else:
            var = 0.0
        sum_a += a
        sum_e += e
        sum_var += var
    if den_or == 0 or sum_var == 0 or n_used == 0:
        return None, None, None, n_used
    or_mh = num_or / den_or
    cmh_stat = (abs(sum_a - sum_e)) ** 2 / sum_var
    p = chi2.sf(cmh_stat, df=1)
    return or_mh, cmh_stat, p, n_used
