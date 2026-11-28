##### BEGIN F08 effect modification: sex, Hispanic and Puerto Rican #####
# Protocol 9.3, 11 and 12.2 analyses 19-20 (ANALYSIS-MAP A19, A20). Gate D15: a group gets its own model only with
# >= 100 confirmed events and >= 21 events in each compared exposure category; otherwise rates only.
# Runs only after the freeze.
import numpy as np
import pandas as pd
import aou_core as C
import _config as K
import analysis as A
import outcomes as O

C.begin("F08", "effect modification: sex (A19), Hispanic / Puerto Rican (A20)")
O.freeze_check()
FR = {q: C.load_state(f"frame_{q}") for q in ("Q1", "Q2")}
MS = K.ANALYSIS["mice_m_sensitivity"]
PAIRS = {"Q1": [("single", "recurrent")], "Q2": [("A", "B"), ("A", "C"), ("A", "D")]}


def joint_fit(F, q, expo, levels, mod, covs):
    """IPTW Cox on F with exposure, modifier (0/1) and product terms, and a joint-category version (RERI).
    Pooled over MS imputations. Returns (product-term p, RERI, lo, hi) for each compared level."""
    sets = O.datasets(F, q, expo, "mice", MS)
    bp, vp, bj, vj = [], [], [], []
    for d in sets:
        d = d.reset_index(drop=True)
        w, _ = O.person_weights(d, expo, levels, covs)
        e = O.expo_matrix(d, expo, levels).iloc[:, 0]
        g = d[mod].astype(float)
        X = pd.DataFrame({"e": e, "g": g, "exg": e * g})
        f = O.cox_fit(d, X, w=w.values, strata=K.ANALYSIS["site_stratum"])
        bp.append(f["beta"]); vp.append(f["vcov"])
        Xj = pd.DataFrame({"e1g0": e * (1 - g), "e0g1": (1 - e) * g, "e1g1": e * g})
        fj = O.cox_fit(d, Xj, w=w.values, strata=K.ANALYSIS["site_stratum"])
        bj.append(fj["beta"]); vj.append(fj["vcov"])
    b, V = A.rubin_vec(bp, vp)
    z = b[2] / np.sqrt(V[2, 2]) if V[2, 2] > 0 else np.nan
    p_int = float(2 * A.stats.norm.sf(abs(z)))
    b2, V2 = A.rubin_vec(bj, vj)
    return p_int, A.reri(b2, V2), float(np.exp(b[2]))


def p_sex():
    """A19: HR in men and in women (PS within each sex), and the exposure x sex product term (Wald)."""
    rows, inter = [], []
    for q in ("Q1", "Q2"):
        F = FR[q].loc[FR[q]["sex"].isin(["Male", "Female"])]
        covs = [c for c in O.adj(q) if c != "sex"]
        for s in ("Female", "Male"):
            rows += O.contrast(F.loc[F["sex"].eq(s)], q, "A19", f"{q} {s.lower()}s", covs=covs)
        for ref, lev in PAIRS[q]:
            G = F.loc[F[O.EXPO[q]].isin([ref, lev])].assign(_male=F["sex"].eq("Male").astype(int))
            ev = G.groupby([O.EXPO[q], "_male"])["event"].sum()
            r = {"question": q, "contrast": f"{O.LEVEL_NAMES[lev]} vs {O.LEVEL_NAMES[ref]}", "modifier": "male vs female"}
            if len(ev) == 4 and ev.min() >= O.MIN_EV:
                r["p_interaction"], (r["RERI"], r["RERI_lo"], r["RERI_hi"]), r["HR_ratio_product"] = \
                    joint_fit(G, q, O.EXPO[q], [ref, lev], "_male", covs + ["sex"])
                r["status"] = "shown"
            else:
                r["status"] = "not shown: a sex x exposure cell has <=20 events"
            inter.append(r)
    O.show_rows("sex", "Sex-specific hazard ratios (A19)", rows, note="PS model fitted within each sex")
    C.show("sex_interaction", "Exposure x sex interaction (A19)", pd.DataFrame(inter), counts=[],
           note=f"product-term Wald p and RERI (delta method), IPTW Cox, MI m = {MS}; exploratory, alpha 0.05, no multiplicity claim")
C.part("1 sex (A19)", p_sex)


def groups(F):
    nhw = F["race_eth"].eq("Non-Hispanic White")
    return {"Hispanic vs non-Hispanic White": (F["hispanic"], nhw),
            "Puerto Rican vs non-Hispanic White": (F["puerto_rican"], nhw & ~F["puerto_rican"])}


def p_hisp():
    """A20 (P§11.2-11.4): group-specific HRs, multiplicative (product term) and additive (RERI) interaction and the
    10-year cumulative-incidence difference by group, only where the D15 gate passes; otherwise descriptive."""
    rows, inter, cifs = [], [], []
    for q in ("Q1", "Q2"):
        F = FR[q]
        for ref, lev in PAIRS[q]:
            G0 = F.loc[F[O.EXPO[q]].isin([ref, lev])]
            for name, (g, r0) in groups(G0).items():
                sub = G0.loc[g | r0].assign(_grp=g.loc[g | r0].astype(int))
                ev = sub.groupby([O.EXPO[q], "_grp"])["event"].sum().reindex(
                    pd.MultiIndex.from_product([[ref, lev], [0, 1]]), fill_value=0)
                grp_ev = int(sub.loc[sub["_grp"].eq(1), "event"].sum())
                gate = grp_ev >= K.DECISIONS["D15_subgroup_min_events"] and ev.min() >= K.DECISIONS["D15_subgroup_min_per_cell"]
                tag = f"{q} {O.LEVEL_NAMES[lev]} vs {O.LEVEL_NAMES[ref]}, {name}"
                ir = {"question": q, "contrast": f"{O.LEVEL_NAMES[lev]} vs {O.LEVEL_NAMES[ref]}", "comparison": name,
                      "D15_gate": "passed" if gate else "descriptive only"}
                if gate:
                    label_g = name.split(" vs ")[0]
                    rows += O.contrast(sub.loc[sub["_grp"].eq(1)], q, "A20", f"{tag}: {label_g}", levels=[ref, lev])
                    rows += O.contrast(sub.loc[sub["_grp"].eq(0)], q, "A20", f"{tag}: non-Hispanic White", levels=[ref, lev])
                    ir["p_interaction"], (ir["RERI"], ir["RERI_lo"], ir["RERI_hi"]), ir["HR_ratio_product"] = \
                        joint_fit(sub, q, O.EXPO[q], [ref, lev], "_grp", O.adj(q))
                    for gv, gl in ((1, label_g), (0, "non-Hispanic White")):
                        c = O.cif(sub.loc[sub["_grp"].eq(gv)], q, method="mice", m=MS, levels=[ref, lev])
                        c = c.loc[c["group"].eq(O.LEVEL_NAMES[lev]) & c["years_after_counting_start"].eq(10.0)]
                        cifs.append(c.assign(comparison=name, subgroup=gl))
                inter.append(ir)
    O.show_rows("hisp", "Hispanic / Puerto Rican group-specific hazard ratios (A20, D15 gate passed only)", rows or
                [{"analysis_id": "A20", "analysis": "no group passed the D15 gate", "status": "descriptive only"}])
    C.show("hisp_interaction", "Hispanic / Puerto Rican: multiplicative and additive interaction (A20)",
           pd.DataFrame(inter), counts=[], note="reference for RERI: non-Hispanic White, unexposed (P§11.3); exploratory")
    if cifs:
        C.show("hisp_cif", "10-year cumulative-incidence difference (exposed - reference) by group (A20)",
               pd.concat(cifs, ignore_index=True), counts=[])
C.part("2 Hispanic / Puerto Rican models (A20)", p_hisp)


def p_hisp_rates():
    """A20 descriptive (P§11.2 a-c): events, person-years and rates per exposure group within Hispanic, Puerto Rican
    (island / mainland) and non-Hispanic White participants. Events are masked with the cross-table check
    (aou_core.leak_hidden) against the arm totals shown in F05.rates."""
    out = []
    for q, levels in (("Q1", ["recurrent", "single", "intermediate"]), ("Q2", ["A", "B", "C", "D"])):
        F = FR[q].loc[FR[q][O.EXPO[q]].isin(levels)].reset_index(drop=True)
        grp = {"Hispanic": F["hispanic"], "Puerto Rican": F["puerto_rican"],
               "Puerto Rican, living in Puerto Rico": F["puerto_rican"] & F["pr_residence"].eq("Puerto Rico"),
               "Puerto Rican, living elsewhere": F["puerto_rican"] & F["pr_residence"].eq("US mainland or other"),
               "Non-Hispanic White": F["race_eth"].eq("Non-Hispanic White")}
        M = {("arm", a): F[O.EXPO[q]].eq(a).values for a in levels}
        for gname, gm in grp.items():
            for a in levels:
                M[(gname, a)] = (gm & F[O.EXPO[q]].eq(a)).values
        hidden, cnt = C.leak_hidden(M, F["event"].values, known=[("arm", a) for a in levels])
        for gname, gm in grp.items():
            for a in levels:
                sel = gm & F[O.EXPO[q]].eq(a)
                py = float(F.loc[sel, "t_years"].sum())
                e = cnt[(gname, a)]
                hid = (gname, a) in hidden
                out.append({"question": q, "group": gname, "exposure": O.LEVEL_NAMES[a],
                            "events": "<=20" if hid else e, "person_years": round(py, 0),
                            "rate_per_1000py": "" if hid or py == 0 else 1000 * e / py})
    C.show("hisp_rates", "RCC rates by exposure within Hispanic, Puerto Rican and non-Hispanic White participants (A20)",
           pd.DataFrame(out), counts=["events"], max_rows=200,
           note="events hidden when 1-20 or computable from other shown cells and the F05 arm totals")
C.part("3 Hispanic / Puerto Rican rates (A20)", p_hisp_rates)
C.end()
##### END F08 #####
