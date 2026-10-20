"""LaTeX tables of the article generated from the pooled replays (data/pooled_*.json).

    python exp70_tables.py            # writes figures/tab_pooled.tex, tab_ablation.tex,
                                      # tab_ablation_v12.tex, tab_baselines_compact.tex,
                                      # tab_baselines_full.tex, tab_transfer.tex,
                                      # tab_heldout.tex, tab_estimate.tex, tab_lmpbase.tex,
                                      # and for Section S15 of the supplement
                                      # tab_versions_summary.tex, tab_versions.tex,
                                      # tab_estimate_v11.tex; prints the totals quoted in the text

Version 6 of the article: counts of 1,000 or more carry a thin separator ({,}), the
ablation tables report the distinct accepted settings and the distinct missed ones,
the thresholds of the baselines are written eta, and Table 7 of the article is the
compact table tab_baselines_compact.tex (the full one goes to the supplement).

Version 7 of the article: the tables of the article carry no columns or rows of earlier
versions of the program.  Those comparisons are in Section S15 of the supplement
(tab_versions_summary.tex, tab_versions.tex, tab_estimate_v11.tex).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "figures")

LABEL = {"L2": "L2", "L3": "L3", "L4": "L4", "L5": "L5", "2D-1": "OpenFOAM 2D-1",
         "2D-2": "OpenFOAM 2D-2", "2D-2h": "OpenFOAM 2D-2", "LMP-20": "LAMMPS 8{,}000", "LMP-40": "LAMMPS 64{,}000",
         "LMP-64": "LAMMPS 262{,}144", "GMX-d": "GROMACS double", "GMX-s": "GROMACS mixed",
         "ABINIT": "ABINIT", "QE": "Quantum ESPRESSO"}
OBS = {"energy": "energy", "dissipation": "dissipation", "gradient": "grad.\\ energy",
       "inhibitor_mass": "inhib.\\ mass", "pattern_sharpness": "sharpness",
       "psi_min": "$\\min\\psi$", "y_min": "$y$ of min.", "omega_lid": "lid vort.",
       "cd": "$c_D$", "dp": "$\\Delta p$", "cl": "$c_L$", "cd_max": "$c_{D,\\max}$",
       "cl_max": "$c_{L,\\max}$", "st": "St", "energy_per_ion": "energy",
       "potential": "energy", "force": "force"}
ORDER = list(LABEL)
PLAIN = {k: v.replace("{,}", ",") for k, v in LABEL.items()}
# Revision 4: which 2D-2 study enters Table 3.  "refined" (default when the rows of the
# Hamilton job J2 are present) shows case 2D-2h and drops the revision-3 case 2D-2;
# "old" shows the revision-3 case.  Both are always in data/pooled_<config>.json.
TABLE3_2D2 = os.environ.get("ERRBUDGET_TABLE3_2D2", "refined")


def _select_2d2(R):
    has_h = any(r["case"] == "2D-2h" for r in R)
    drop = "2D-2" if (has_h and TABLE3_2D2 == "refined") else "2D-2h"
    return [r for r in R if r["case"] != drop]


def load(cfg, select=True):
    p = os.path.join(D, f"pooled_{cfg}.json")
    if not os.path.exists(p):
        return None
    R = json.load(open(p))
    if select:
        R = _select_2d2(R)
    return sorted(R, key=lambda r: (ORDER.index(r["case"]), list(OBS).index(r["ob"]), -r["target"]))


def tex_num(x):
    if x is None or not np.isfinite(x):
        return "---"
    return f"{x:.2f}" if x < 10 else (f"{x:.1f}" if x < 100 else f"{x:.0f}")


def tex_int(n):
    """An integer with a thin separator for thousands, as in the article (17{,}610)."""
    return f"{int(n):,}".replace(",", "{,}")


def tex_big(x):
    """tex_num with a thousands separator for values of 1,000 or more."""
    s = tex_num(x)
    return tex_int(round(x)) if s != "---" and x >= 1000 else s


def tex_tg(t):
    e = int(np.floor(np.log10(t) + 1e-9))
    m = round(t / 10**e)
    if m == 10:
        e, m = e + 1, 1
    return f"$10^{{{e}}}$" if m == 1 else f"${m}{{\\times}}10^{{{e}}}$"


def totals(R):
    t = lambda k: sum(r.get(k, 0) for r in R)
    pm = [r["pen_med"] for r in R if r["accepted"] and np.isfinite(r["pen_med"])]
    return dict(cells=len(R), n=t("n"), sf=t("start_failed"), acc=t("accepted"), miss=t("missed"),
                dr=t("declined_reach"), du=t("declined_unreach"), unconv=t("unconv"),
                distinct=t("distinct"), distinct_missed=t("distinct_missed"),
                cells_miss=sum(1 for r in R if r["missed"]),
                pen=float(np.median(pm)) if pm else np.nan,
                geo=float(np.exp(np.mean(np.log(pm)))) if pm else np.nan,
                runs=float(np.median([r["runs_med"] for r in R if np.isfinite(r["runs_med"])])))


def pooled_table(R):
    lines = []
    last = None
    for r in R:
        lab = LABEL[r["case"]] if r["case"] != last else ""
        if r["case"] != last and last is not None:
            lines.append("\\addlinespace[1pt]")
        last = r["case"]
        lines.append(f"{lab} & {OBS[r['ob']]} & {tex_tg(r['target'])} & {tex_int(r['n'])} & {tex_int(r['start_failed'])} & "
                     f"{tex_int(r['accepted'])} & {r.get('distinct', '---')} & {tex_int(r['missed'])} & "
                     f"{tex_int(r['declined_reach'])} & {tex_int(r['declined_unreach'])} & "
                     f"{tex_num(r['pen_med'])} & {tex_num(r['pen_max'])} & {r['runs_med']:.0f} \\\\")
    T = totals(R)
    lines.append("\\midrule")
    lines.append(f"\\multicolumn{{3}}{{l}}{{All {T['cells']} cells}} & {tex_int(T['n'])} & {tex_int(T['sf'])} & "
                 f"{tex_int(T['acc'])} & {tex_int(T['distinct'])} & {tex_int(T['miss'])} & {tex_int(T['dr'])} & "
                 f"{tex_int(T['du'])} & {tex_num(T['pen'])} & & {T['runs']:.0f} \\\\")
    return "\n".join(lines)


def versions_table(cfgs):
    """Section S15: accepted/missed in every cell of Table 3 for each version of the program."""
    data = [(c, {(r["case"], r["ob"], r["target"]): r for r in (load(c) or [])}) for c, _ in cfgs]
    base = load(cfgs[0][0]) or []
    lines, last = [], None
    for r in base:
        lab = LABEL[r["case"]] if r["case"] != last else ""
        if r["case"] != last and last is not None:
            lines.append("\\addlinespace[1pt]")
        last = r["case"]
        k = (r["case"], r["ob"], r["target"])
        cells = []
        for _, D_ in data:
            o = D_.get(k)
            cells.append(f"{tex_int(o['accepted'])}/{tex_int(o['missed'])}" if o else "---")
        lines.append(f"{lab} & {OBS[r['ob']]} & {tex_tg(r['target'])} & {tex_int(r['n'])} & "
                     + " & ".join(cells) + " \\\\")
    tots = []
    for c, D_ in data:
        T = totals(list(D_.values()))
        tots.append(f"{tex_int(T['acc'])}/{tex_int(T['miss'])}")
    lines.append("\\midrule")
    lines.append(f"\\multicolumn{{3}}{{l}}{{All {len(base)} cells}} & {tex_int(totals(base)['n'])} & "
                 + " & ".join(tots) + " \\\\")
    return "\n".join(lines)


def ablation_table(cfgs, base_cfg="v1.2.1"):
    rows = []
    base = load(base_cfg)
    keys = {(r["case"], r["ob"], r["target"]) for r in base} if base else set()
    for name, lab in cfgs:
        R = load(name)
        if R is None:
            continue
        R = [r for r in R if (r["case"], r["ob"], r["target"]) in keys]
        T = totals(R)
        cells = sorted({r["case"] for r in R if r["missed"]}, key=ORDER.index)
        cells = list(dict.fromkeys(LABEL[c] for c in cells))
        miss = tex_int(T["miss"]) + (f" ({T['distinct_missed']})" if T["miss"] else "")
        rows.append(f"{lab} & {tex_int(T['acc'])} & {tex_int(T['distinct'])} & {miss} & {T['cells_miss']} & "
                    f"{tex_int(T['dr'])} & {tex_num(T['pen'])} & {tex_num(T['geo'])} & {T['runs']:.0f} & "
                    f"{', '.join(cells) if cells else '---'} \\\\")
    return "\n".join(rows)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    R = load("v1.2.1")
    open(os.path.join(OUT, "tab_pooled.tex"), "w").write(pooled_table(R))
    for name in ["v1.2.1", "v1.2", "v1.1"]:
        X = load(name)
        if X:
            print(name, totals(X))
    # Revision 4: every component is removed from version 1.2.1, the configuration of
    # Table 3 (revision 3 removed them from the frozen version 1.2; that table is kept
    # as tab_ablation_v12.tex for the response letter).
    # Version 7: Table 5 of the article removes components from the released program and
    # shows the frozen program as its joint-check threshold; the earlier versions are in
    # Section S15 of the supplement (tab_versions_summary.tex).
    abl = [("v1.2.1", "Default configuration"),
           ("v1.2", "joint check only above $3\\rho$ (frozen program)"),
           ("v1.2.1-nojoint", "without the joint check"),
           ("v1.2.1-nolocal", "without the local-order check"),
           ("v1.2.1-credit", "crediting cancellation"),
           ("v1.2.1-constant", "constant safety factor"),
           ("v1.2.1-nores", "resolution without the $0.1\\,\\varepsilon$ term"),
           ("v1.2.1-accept-unconverged", "accepting unconverged runs"),
           ("v1.2.1-nested-as-disc", "nested solver declared discretization"),
           ("v1.2.1-screened", "with the preliminary sweep")]
    open(os.path.join(OUT, "tab_ablation.tex"), "w").write(ablation_table(abl))
    vers = [("v1.2.1", "1.2.1 (released)"),
            ("v1.2", "1.2 (frozen)"),
            ("v1.1+gated", "1.1, gated factor"),
            ("v1.1", "1.1")]
    open(os.path.join(OUT, "tab_versions_summary.tex"), "w").write(ablation_table(vers))
    open(os.path.join(OUT, "tab_versions.tex"), "w").write(versions_table(vers))
    abl12 = [(n.replace("v1.2.1-", "v1.2-"), l) for n, l in abl if n.startswith("v1.2.1-")]
    open(os.path.join(OUT, "tab_ablation_v12.tex"), "w").write(
        ablation_table([("v1.2", "Version 1.2 (frozen)")] + abl12, base_cfg="v1.2"))
    print(open(os.path.join(OUT, "tab_ablation.tex")).read())


def baseline_cells():
    p = os.path.join(D, "baselines_v6.json")
    if os.path.exists(p):
        res = json.load(open(p))["cells"]
    else:
        p = os.path.join(D, "baselines_pooled.json")
        if not os.path.exists(p):
            return None, []
        res = json.load(open(p))
    has_h = any(k.split("|")[0] == "2D-2h" for k in res)
    two = (["2D-2h"] if (has_h and TABLE3_2D2 == "refined") else ["2D-2"])
    groups = [("L2--L4", ["L2", "L3", "L4"]), ("L5", ["L5"]), ("OpenFOAM 2D-1", ["2D-1"]),
              ("OpenFOAM 2D-2", two), ("LAMMPS", ["LMP-20", "LMP-40", "LMP-64"]),
              ("GROMACS double", ["GMX-d"]), ("GROMACS mixed", ["GMX-s"]), ("ABINIT", ["ABINIT"]),
              ("Quantum ESPRESSO", ["QE"])]
    return res, groups


def baselines_compact_table():
    """Table 7 of the article: one row per group.  Program: cells, starts, declines, penalty,
    runs and K; misses of sequential convergence (S) and coordinate search (C) at the three
    thresholds and of reference acceptance; penalty and K of C at eta = eps/3 and eps/10.
    Group values are medians over the cells of the per-cell medians, counts are sums."""
    res, groups = baseline_cells()
    if res is None:
        return ""
    med = lambda R, f: float(np.median([r[f] for r in R if np.isfinite(r[f])])) if any(
        np.isfinite(r[f]) for r in R) else np.nan
    out = []
    for g, ks in groups:
        rows = [v for k, v in res.items() if k.split("|")[0] in ks]
        if not rows:
            continue
        P = [v["program"] for v in rows]
        starts = sum(r["reached"] + r["missed"] + r["declined"] for r in P)
        cells = [g, str(len(rows)), tex_int(starts), tex_int(sum(r["declined"] for r in P)),
                 tex_num(med(P, "pen")), f"{np.median([r['runs'] for r in P]):.0f}", tex_big(med(P, "K"))]
        for key in ("seq-1", "seq-3", "seq-10", "coord-1", "coord-3", "coord-10", "ref"):
            cells.append(tex_int(sum(v[key]["missed"] for v in rows)))
        for key in ("coord-3", "coord-10"):
            R = [v[key] for v in rows]
            cells += [tex_num(med(R, "pen")), f"{med(R, 'K'):.1f}"]
        out.append(" & ".join(cells) + " \\\\")
    return "\n".join(out)


def baselines_table():
    res, groups = baseline_cells()
    if res is None:
        return ""
    names = [("program", "This program"), ("seq-1", "Sequential, $\\eta=\\varepsilon$"),
             ("seq-3", "Sequential, $\\eta=\\varepsilon/3$"), ("seq-10", "Sequential, $\\eta=\\varepsilon/10$"),
             ("coord-1", "Coordinate, $\\eta=\\varepsilon$"), ("coord-3", "Coordinate, $\\eta=\\varepsilon/3$"),
             ("coord-10", "Coordinate, $\\eta=\\varepsilon/10$"), ("ref", "Reference acceptance")]
    out = []
    for g, ks in groups:
        rows = {k: v for k, v in res.items() if k.split("|")[0] in ks}
        if not rows:
            continue
        first = True
        for key, lab in names:
            R = [v[key] for v in rows.values() if key in v]
            if not R:
                continue
            pen = [r["pen"] for r in R if np.isfinite(r["pen"])]
            Ks = [r["K"] for r in R if np.isfinite(r["K"])]
            dec = sum(r["declined"] for r in R)
            out.append(f"{g + f' ({len(rows)} cells)' if first else ''} & {lab} & {sum(r['reached'] for r in R)} & "
                       f"{sum(r['missed'] for r in R)} & {dec if key == 'program' else '---'} & "
                       f"{tex_num(np.median(pen)) if pen else '---'} & {np.median([r['runs'] for r in R]):.0f} & "
                       f"{tex_num(np.median(Ks)) if Ks else '---'} \\\\")
            first = False
        out.append("\\addlinespace[2pt]")
    return "\n".join(out[:-1])


open(os.path.join(OUT, "tab_baselines_full.tex"), "w").write(baselines_table())
open(os.path.join(OUT, "tab_baselines_compact.tex"), "w").write(baselines_compact_table())


def transfer_table():
    p = os.path.join(D, "transfer_floor.json")
    if not os.path.exists(p):
        return ""
    rows = json.load(open(p))
    cls = {"discretisation": "d", "nested-solver": "n", "cost-only": "c"}
    lab = {"L2": "L2", "L3": "L3", "L4": "L4", "L5": "L5", "2D-1": "2D-1", "2D-2": "2D-2",
           "QE": "QE", "ABINIT": "ABINIT"}
    out = ["\\begin{tabular}{lllrrrrrrrrrr}", "\\toprule",
           "& & & \\multicolumn{5}{c}{Without the refinements} & \\multicolumn{5}{c}{With the refinements} \\\\",
           "\\cmidrule(lr){4-8}\\cmidrule(lr){9-13}",
           "Problem & Parameter & & $f{=}0.01$ & $0.05$ & $0.1$ & $0.25$ & $0.5$ & $0.01$ & $0.05$ & $0.1$ & $0.25$ & $0.5$ \\\\",
           "\\midrule"]
    last = None
    keyed = {}
    for r in rows:
        keyed.setdefault((r["case"], r["param"], r["kind"]), {})[r["floor"]] = r["rates"]
    for (case, par, kind), v in keyed.items():
        if last is not None and case != last:
            out.append("\\addlinespace[1pt]")
        a = v.get(False, [np.nan] * 5); b = v.get(True, [np.nan] * 5)
        pname = par.replace("_", "\\_")
        out.append(f"{lab[case] if case != last else ''} & {pname} & {cls[kind]} & "
                   + " & ".join(f"{x:.0f}\\%" for x in a) + " & " + " & ".join(f"{x:.0f}\\%" for x in b) + " \\\\")
        last = case
    out += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(out)


open(os.path.join(OUT, "tab_transfer.tex"), "w").write(transfer_table())


def heldout_table():
    p1 = os.path.join(D, "heldout_replay.json"); p2 = os.path.join(D, "heldout_replay_current.json")
    if not os.path.exists(p1):
        return ""
    R = json.load(open(p1)); Rc = json.load(open(p2)) if os.path.exists(p2) else []
    get = lambda rows, cfg: {(r["case"], r["ob"], r["target"]): r for r in rows if r["config"] == cfg}
    v12, v121 = get(R, "v1.2"), get(Rc, "v1.2.1")
    obs = {"mean": "mean", "flux": "flux", "u_min": "$u_{\\min}$", "v_max": "$v_{\\max}$",
           "v_min": "$v_{\\min}$", "ek": "$E_k$", "eps": "$\\epsilon$", "umax": "$u_{\\max}$"}
    out, last = [], None
    for k, r in v12.items():
        case = k[0]
        if last is not None and case != last:
            out.append("\\addlinespace[1pt]")
        c = v121.get(k)
        out.append(f"{case if case != last else ''} & {obs.get(k[1], k[1])} & {tex_tg(k[2])} & {r['n']} & "
                   f"{r['start_failed']} & {r['accepted']} & {r.get('distinct', '---')} & {r['missed']} & {r['declined_reach'] + r['declined_unreach']} & "
                   f"{tex_num(r['pen_med'])} & {r['runs_med']:.0f} & "
                   f"{(str(c['accepted']) + '/' + str(c['missed'])) if c else '---'} \\\\")
        last = case
    return "\n".join(out)


open(os.path.join(OUT, "tab_heldout.tex"), "w").write(heldout_table())


def estimate_table(version="v12"):
    """Table A.1 of the article (version="v12": the released program, whose results the
    file stores under the key v12) or the version-1.1 columns of Section S15 ("v11")."""
    p = os.path.join(D, "estimate_v12.json")
    if not os.path.exists(p):
        return ""
    R = json.load(open(p))
    head = "Final estimate" if version == "v12" else "Version 1.1"
    out = ["\\begin{tabular}{llrrrrr}", "\\toprule",
           f"& & & Raw & \\multicolumn{{3}}{{c}}{{{head}}} \\\\",
           "\\cmidrule(lr){5-7}",
           "Problem & Observable & Settings & opt. & 10th pct. & Median & Opt. \\\\",
           "\\midrule"]
    last = None
    for r in R:
        lab = LABEL.get(r["case"], r["case"]) if r["case"] != last else ""
        last = r["case"]
        a = r[version]
        out.append(f"{lab} & {OBS.get(r['ob'], r['ob'])} & {r['n']} & {r['raw']:.0f}\\% & {a[0]:.2f} & {a[1]:.2f} & "
                   f"{a[2]:.1f}\\% \\\\")
    out += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(out)


open(os.path.join(OUT, "tab_estimate.tex"), "w").write(estimate_table("v12"))
open(os.path.join(OUT, "tab_estimate_v11.tex"), "w").write(estimate_table("v11"))


def lmp_table():
    p = os.path.join(D, "baselines_lmp_gmx.json")
    if not os.path.exists(p):
        return ""
    R = json.load(open(p))
    if isinstance(R, dict) and "cells" in R:          # revision 4 format of exp72: {config, cells}
        R = R["cells"]
    lab = {"LMP-20": "LAMMPS 8,000", "LMP-40": "LAMMPS 64,000", "LMP-64": "LAMMPS 262,144",
           "GMX-d": "GROMACS double"}
    procs = [("program", "This program"), ("seq-1", "Sequential, $\\eta=\\varepsilon$"),
             ("seq-10", "Sequential, $\\eta=\\varepsilon/10$"),
             ("coord-10", "Coordinate, $\\eta=\\varepsilon/10$"), ("ref", "Reference acceptance")]
    out = ["\\begin{tabular}{ll" + "r" * len(procs) + "}", "\\toprule",
           "Case & Target & " + " & ".join(p[1] for p in procs) + " \\\\", "\\midrule"]
    for k, v in R.items():
        c, ob, tg = k.split("|")
        cells = []
        for p, _ in procs:
            x = v[p]; n = x["reached"] + x["missed"] + x["declined"]
            pen = "---" if not np.isfinite(x["pen"]) else f"{x['pen']:.2f}"
            cells.append(f"{x['reached']}/{n} ({pen})")
        out.append(f"{lab[c]} & {tex_tg(float(tg))} & " + " & ".join(cells) + " \\\\")
    out += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(out)


open(os.path.join(OUT, "tab_lmpbase.tex"), "w").write(lmp_table())
