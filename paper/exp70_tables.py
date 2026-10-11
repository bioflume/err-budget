"""LaTeX tables of the article generated from the pooled replays (data/pooled_*.json).

    python exp70_tables.py            # writes figures/tab_pooled.tex, tab_ablation.tex,
                                      # tab_versions.tex and prints the totals quoted in the text
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "figures")

LABEL = {"L2": "L2", "L3": "L3", "L4": "L4", "L5": "L5", "2D-1": "OpenFOAM 2D-1",
         "2D-2": "OpenFOAM 2D-2", "LMP-20": "LAMMPS 8{,}000", "LMP-40": "LAMMPS 64{,}000",
         "LMP-64": "LAMMPS 262{,}144", "GMX-d": "GROMACS double", "GMX-s": "GROMACS mixed",
         "ABINIT": "ABINIT", "QE": "Quantum ESPRESSO"}
OBS = {"energy": "energy", "dissipation": "dissipation", "gradient": "grad.\\ energy",
       "inhibitor_mass": "inhib.\\ mass", "pattern_sharpness": "sharpness",
       "psi_min": "$\\min\\psi$", "y_min": "$y$ of min.", "omega_lid": "lid vort.",
       "cd": "$c_D$", "dp": "$\\Delta p$", "cl": "$c_L$", "cd_max": "$c_{D,\\max}$",
       "cl_max": "$c_{L,\\max}$", "st": "St", "energy_per_ion": "energy",
       "potential": "energy", "force": "force"}
ORDER = list(LABEL)


def load(cfg):
    p = os.path.join(D, f"pooled_{cfg}.json")
    if not os.path.exists(p):
        return None
    R = json.load(open(p))
    return sorted(R, key=lambda r: (ORDER.index(r["case"]), list(OBS).index(r["ob"]), -r["target"]))


def tex_num(x):
    if x is None or not np.isfinite(x):
        return "---"
    return f"{x:.2f}" if x < 10 else (f"{x:.1f}" if x < 100 else f"{x:.0f}")


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
                cells_miss=sum(1 for r in R if r["missed"]),
                pen=float(np.median(pm)) if pm else np.nan,
                geo=float(np.exp(np.mean(np.log(pm)))) if pm else np.nan,
                runs=float(np.median([r["runs_med"] for r in R if np.isfinite(r["runs_med"])])))


def pooled_table(R, R11):
    by11 = {(r["case"], r["ob"], r["target"]): r for r in (R11 or [])}
    lines = []
    last = None
    for r in R:
        lab = LABEL[r["case"]] if r["case"] != last else ""
        if r["case"] != last and last is not None:
            lines.append("\\addlinespace[1pt]")
        last = r["case"]
        o = by11.get((r["case"], r["ob"], r["target"]))
        v11 = f"{o['accepted']}/{o['missed']}" if o else "---"
        lines.append(f"{lab} & {OBS[r['ob']]} & {tex_tg(r['target'])} & {r['n']} & {r['start_failed']} & "
                     f"{r['accepted']} & {r['missed']} & {r['declined_reach']} & {r['declined_unreach']} & "
                     f"{tex_num(r['pen_med'])} & {tex_num(r['pen_max'])} & {r['runs_med']:.0f} & {v11} \\\\")
    T = totals(R)
    T11 = totals(R11) if R11 else None
    lines.append("\\midrule")
    lines.append(f"\\multicolumn{{3}}{{l}}{{All {T['cells']} cells}} & {T['n']} & {T['sf']} & {T['acc']} & "
                 f"{T['miss']} & {T['dr']} & {T['du']} & {tex_num(T['pen'])} & & {T['runs']:.0f} & "
                 + (f"{T11['acc']}/{T11['miss']}" if T11 else "---") + " \\\\")
    return "\n".join(lines)


def ablation_table(cfgs):
    rows = []
    base = load("v1.2")
    keys = {(r["case"], r["ob"], r["target"]) for r in base} if base else set()
    for name, lab in cfgs:
        R = load(name)
        if R is None:
            continue
        R = [r for r in R if (r["case"], r["ob"], r["target"]) in keys]
        T = totals(R)
        cells = sorted({r["case"] for r in R if r["missed"]}, key=ORDER.index)
        rows.append(f"{lab} & {T['acc']} & {T['miss']} & {T['cells_miss']} & {T['dr']} & "
                    f"{tex_num(T['pen'])} & {tex_num(T['geo'])} & {T['runs']:.0f} & "
                    f"{', '.join(cells) if cells else '---'} \\\\")
    return "\n".join(rows)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    R, R11 = load("v1.2.1"), load("v1.1")
    open(os.path.join(OUT, "tab_pooled.tex"), "w").write(pooled_table(R, R11))
    for name in ["v1.2.1", "v1.2", "v1.1"]:
        X = load(name)
        if X:
            print(name, totals(X))
    abl = [("v1.2.1", "Version 1.2.1 (default)"),
           ("v1.2", "Version 1.2 (frozen)"),
           ("v1.2-nojoint", "without the joint check"),
           ("v1.2-nolocal", "without the local-order check"),
           ("v1.2-credit", "crediting cancellation"),
           ("v1.2-constant", "constant safety factor"),
           ("v1.2-nores", "without the resolution floor"),
           ("v1.2-accept-unconverged", "accepting unconverged runs"),
           ("v1.2-nested-as-disc", "nested solver declared discretisation"),
           ("v1.2-screened", "with the preliminary sweep"),
           ("v1.1+gated", "version 1.1 with gated factor"),
           ("v1.1", "version 1.1")]
    open(os.path.join(OUT, "tab_ablation.tex"), "w").write(ablation_table(abl))
    print(open(os.path.join(OUT, "tab_ablation.tex")).read())


def baselines_table():
    p = os.path.join(D, "baselines_pooled.json")
    if not os.path.exists(p):
        return ""
    res = json.load(open(p))
    groups = [("L2--L4", ["L2", "L3", "L4"]), ("L5", ["L5"]), ("OpenFOAM 2D-1", ["2D-1"]),
              ("OpenFOAM 2D-2", ["2D-2"]), ("Quantum ESPRESSO", ["QE"])]
    names = [("program", "This program"), ("seq-1", "Sequential, $\\delta=\\varepsilon$"),
             ("seq-3", "Sequential, $\\delta=\\varepsilon/3$"), ("seq-10", "Sequential, $\\delta=\\varepsilon/10$"),
             ("coord-10", "Coordinate, $\\delta=\\varepsilon/10$"), ("ref", "Reference acceptance")]
    out = []
    for g, ks in groups:
        rows = {k: v for k, v in res.items() if k.split("|")[0] in ks}
        if not rows:
            continue
        first = True
        for key, lab in names:
            R = [v[key] for v in rows.values()]
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


open(os.path.join(OUT, "tab_baselines.tex"), "w").write(baselines_table())


def transfer_table():
    p = os.path.join(D, "transfer_floor.json")
    if not os.path.exists(p):
        return ""
    rows = json.load(open(p))
    cls = {"discretisation": "d", "nested-solver": "n", "cost-only": "c"}
    lab = {"L2": "L2", "L3": "L3", "L4": "L4", "L5": "L5", "2D-1": "2D-1", "2D-2": "2D-2",
           "QE": "QE", "ABINIT": "ABINIT"}
    out = ["\\begin{tabular}{lllrrrrrrrrrr}", "\\toprule",
           "& & & \\multicolumn{5}{c}{Without floor} & \\multicolumn{5}{c}{With floor} \\\\",
           "\\cmidrule(lr){4-8}\\cmidrule(lr){9-13}",
           "Problem & Parameter & & $\\phi{=}0.01$ & $0.05$ & $0.1$ & $0.25$ & $0.5$ & $0.01$ & $0.05$ & $0.1$ & $0.25$ & $0.5$ \\\\",
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
    v12, v11, v121 = get(R, "v1.2"), get(R, "v1.1"), get(Rc, "v1.2.1")
    obs = {"mean": "mean", "flux": "flux", "u_min": "$u_{\\min}$", "v_max": "$v_{\\max}$",
           "v_min": "$v_{\\min}$", "ek": "$E_k$", "eps": "$\\epsilon$", "umax": "$u_{\\max}$"}
    out, last = [], None
    for k, r in v12.items():
        case = k[0]
        if last is not None and case != last:
            out.append("\\addlinespace[1pt]")
        o, c = v11.get(k), v121.get(k)
        out.append(f"{case if case != last else ''} & {obs.get(k[1], k[1])} & {tex_tg(k[2])} & {r['n']} & "
                   f"{r['start_failed']} & {r['accepted']} & {r['missed']} & {r['declined_reach'] + r['declined_unreach']} & "
                   f"{tex_num(r['pen_med'])} & {r['runs_med']:.0f} & "
                   f"{(str(o['accepted']) + '/' + str(o['missed'])) if o else '---'} & "
                   f"{(str(c['accepted']) + '/' + str(c['missed'])) if c else '---'} \\\\")
        last = case
    return "\n".join(out)


open(os.path.join(OUT, "tab_heldout.tex"), "w").write(heldout_table())


def estimate_table():
    p = os.path.join(D, "estimate_v12.json")
    if not os.path.exists(p):
        return ""
    R = json.load(open(p))
    out = ["\\begin{tabular}{llrrrrrrrr}", "\\toprule",
           "& & & Raw & \\multicolumn{3}{c}{Version 1.1} & \\multicolumn{3}{c}{Version 1.2.1} \\\\",
           "\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
           "Problem & Observable & Settings & opt. & 10th pct. & Median & Opt. & 10th pct. & Median & Opt. \\\\",
           "\\midrule"]
    last = None
    for r in R:
        lab = LABEL.get(r["case"], r["case"]) if r["case"] != last else ""
        last = r["case"]
        a, b = r["v11"], r["v12"]
        out.append(f"{lab} & {OBS.get(r['ob'], r['ob'])} & {r['n']} & {r['raw']:.0f}\\% & {a[0]:.2f} & {a[1]:.2f} & "
                   f"{a[2]:.1f}\\% & {b[0]:.2f} & {b[1]:.2f} & {b[2]:.1f}\\% \\\\")
    out += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(out)


open(os.path.join(OUT, "tab_estimate.tex"), "w").write(estimate_table())


def lmp_table():
    p = os.path.join(D, "baselines_lmp_gmx.json")
    if not os.path.exists(p):
        return ""
    R = json.load(open(p))
    lab = {"LMP-20": "LAMMPS 8,000", "LMP-40": "LAMMPS 64,000", "LMP-64": "LAMMPS 262,144",
           "GMX-d": "GROMACS double"}
    procs = [("program", "This program"), ("seq-1", "Sequential, $\\delta=\\varepsilon$"),
             ("seq-10", "Sequential, $\\delta=\\varepsilon/10$"),
             ("coord-10", "Coordinate, $\\delta=\\varepsilon/10$"), ("ref", "Reference acceptance")]
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
