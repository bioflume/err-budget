"""The LaTeX tables of v12, from the stored replays in data/.

    python exp70_tables.py     # writes figures/tab_*.tex and prints the totals quoted in the text

Article: Table 4 (tab_pooled_compact.tex), Table 5 (tab_preregistered.tex).
Supplement: Tables S-3.1 (tab_pooled.tex), S-3.2 (tab_sup_ablation.tex), S-3.3
(tab_estimate.tex), S-4.1 (tab_transfer.tex), S-4.2 (tab_baselines_full.tex) and S-4.3
(tab_lmpbase.tex).  Table S-3.5 (tab_misdeclare.tex) is written by exp65_misdeclare_v12.py.

Group values are medians over the test cases of the per-case medians; counts are sums.
2D-2 is the refined study (case 2D-2h); the coarser study (case 2D-2) is left out.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "figures")

LABEL = {"L2": "L2", "L3": "L3", "L4": "L4", "L5": "L5", "2D-1": "OpenFOAM 2D-1",
         "2D-2h": "OpenFOAM 2D-2", "LMP-20": "LAMMPS 8{,}000", "LMP-40": "LAMMPS 64{,}000",
         "LMP-64": "LAMMPS 262{,}144", "GMX-d": "GROMACS double", "GMX-s": "GROMACS mixed",
         "ABINIT": "ABINIT", "QE": "Quantum ESPRESSO"}
OBS = {"energy": "energy", "dissipation": "dissipation", "gradient": "grad.\\ energy",
       "inhibitor_mass": "inhib.\\ mass", "pattern_sharpness": "sharpness",
       "psi_min": "$\\min\\psi$", "y_min": "$y$ of min.", "omega_lid": "lid vort.",
       "cd": "$c_D$", "dp": "$\\Delta p$", "cl": "$c_L$", "cd_max": "$c_{D,\\max}$",
       "cl_max": "$c_{L,\\max}$", "st": "St", "energy_per_ion": "energy",
       "potential": "energy", "force": "force"}
GROUPS = [("L2--L4", ["L2", "L3", "L4"]), ("L5", ["L5"]), ("OpenFOAM 2D-1", ["2D-1"]),
          ("OpenFOAM 2D-2", ["2D-2h"]), ("LAMMPS", ["LMP-20", "LMP-40", "LMP-64"]),
          ("GROMACS double", ["GMX-d"]), ("GROMACS mixed", ["GMX-s"]), ("ABINIT", ["ABINIT"]),
          ("Quantum ESPRESSO", ["QE"])]
PROCEDURES = [("program", "This program"), ("seq-1", "S, $\\eta=\\varepsilon$"),
              ("seq-3", "S, $\\eta=\\varepsilon/3$"), ("seq-10", "S, $\\eta=\\varepsilon/10$"),
              ("coord-1", "SJ, $\\eta=\\varepsilon$"), ("coord-3", "SJ, $\\eta=\\varepsilon/3$"),
              ("coord-10", "SJ, $\\eta=\\varepsilon/10$"), ("ref", "Ref.")]
ABLATION = [("v1.2.1", "Default configuration"),
            ("v1.2.1-nojoint", "\\quad without the interaction check"),
            ("v1.2.1-nolocal", "\\quad without the observed-order check"),
            ("v1.2.1-credit", "\\quad signed sum (cancellation allowed)"),
            ("v1.2.1-constant", "\\quad fixed safety factor 1.25"),
            ("v1.2.1-nores", "\\quad detection threshold $\\rho=\\sigma$"),
            ("v1.2.1-accept-unconverged", "\\quad accepting unconverged runs"),
            ("v1.2.1-nested-as-disc", "\\quad iterative-solver parameters declared discretization"),
            ("v1.2.1-screened", "\\quad with the screening sweep")]


def load_json(name):
    return json.load(open(os.path.join(D, name)))


def load(cfg):
    """Rows of data/pooled_<cfg>.json in the order of the tables, without the coarser 2D-2."""
    R = [r for r in load_json(f"pooled_{cfg}.json") if r["case"] in LABEL]
    return sorted(R, key=lambda r: (list(LABEL).index(r["case"]), list(OBS).index(r["ob"]),
                                    -r["target"]))


def tex_num(x):
    if x is None or not np.isfinite(x):
        return "---"
    return f"{x:.2f}" if x < 10 else (f"{x:.1f}" if x < 100 else f"{x:.0f}")


def tex_int(n):
    """An integer with a thin separator for thousands (17{,}610)."""
    return f"{int(n):,}".replace(",", "{,}")


def tex_tg(t):
    e = int(np.floor(np.log10(t) + 1e-9))
    m = round(t / 10**e)
    if m == 10:
        e, m = e + 1, 1
    return f"$10^{{{e}}}$" if m == 1 else f"${m}{{\\times}}10^{{{e}}}$"


def median(xs):
    xs = [x for x in xs if np.isfinite(x)]
    return float(np.median(xs)) if xs else np.nan


def totals(R):
    t = lambda k: sum(r.get(k, 0) for r in R)
    pm = [r["pen_med"] for r in R if r["accepted"] and np.isfinite(r["pen_med"])]
    return dict(cases=len(R), n=t("n"), sf=t("start_failed"), acc=t("accepted"), miss=t("missed"),
                dr=t("declined_reach"), du=t("declined_unreach"), unconv=t("unconv"),
                distinct=t("distinct"), distinct_missed=t("distinct_missed"),
                cases_miss=sum(1 for r in R if r["missed"]),
                pen=float(np.median(pm)) if pm else np.nan,
                geo=float(np.exp(np.mean(np.log(pm)))) if pm else np.nan,
                runs=median([r["runs_med"] for r in R]))


def pooled_table(R):
    """Table S-3.1: every development test case."""
    lines, last = [], None
    for r in R:
        if r["case"] != last and last is not None:
            lines.append("\\addlinespace[1pt]")
        lab = LABEL[r["case"]] if r["case"] != last else ""
        last = r["case"]
        lines.append(f"{lab} & {OBS[r['ob']]} & {tex_tg(r['target'])} & {tex_int(r['n'])} & "
                     f"{tex_int(r['start_failed'])} & {tex_int(r['accepted'])} & {r['distinct']} & "
                     f"{tex_int(r['missed'])} & {tex_int(r['declined_reach'])} & "
                     f"{tex_int(r['declined_unreach'])} & {tex_num(r['pen_med'])} & "
                     f"{tex_num(r['pen_max'])} & {r['runs_med']:.0f} \\\\")
    T = totals(R)
    lines += ["\\midrule",
              f"\\multicolumn{{3}}{{l}}{{All {T['cases']} test cases}} & {tex_int(T['n'])} & "
              f"{tex_int(T['sf'])} & {tex_int(T['acc'])} & {tex_int(T['distinct'])} & "
              f"{tex_int(T['miss'])} & {tex_int(T['dr'])} & {tex_int(T['du'])} & "
              f"{tex_num(T['pen'])} & & {T['runs']:.0f} \\\\"]
    return "\n".join(lines)


def compact_table(R, base):
    """Table 4: the development test cases by group.  Cost ratio and runs per search are
    those of the program from converged starting settings (data/baselines_v6.json)."""
    lines = []
    for g, keys in GROUPS:
        rows = [r for r in R if r["case"] in keys]
        prog = [v["program"] for k, v in base.items() if k.split("|")[0] in keys]
        T = totals(rows)
        lines.append(f"{g} & {T['cases']} & {tex_int(T['n'])} & {tex_int(T['acc'])} & "
                     f"{T['distinct']} & {T['miss']} & {tex_int(T['dr'] + T['du'])} & "
                     f"{tex_num(median([p['pen'] for p in prog]))} & "
                     f"{np.median([p['runs'] for p in prog]):.0f} \\\\")
    T = totals(R)
    prog = [v["program"] for k, v in base.items() if k.split("|")[0] in LABEL]
    lines += ["\\midrule",
              f"All groups & {T['cases']} & {tex_int(T['n'])} & {tex_int(T['acc'])} & "
              f"{T['distinct']} & {T['miss']} & {tex_int(T['dr'] + T['du'])} & "
              f"{tex_num(T['pen'])} & {np.median([p['runs'] for p in prog]):.0f} \\\\"]
    return "\n".join(lines) + "\n"


def preregistered_table():
    """Table 5: the searches of the frozen program on H1-H3 and the released program's
    accepted and missed searches."""
    frozen = {(r["case"], r["ob"], r["target"]): r
              for r in load_json("heldout_replay.json") if r["config"] == "v1.2"}
    released = {(r["case"], r["ob"], r["target"]): r
                for r in load_json("heldout_replay_current.json") if r["config"] == "v1.2.1"}
    obs = {"mean": "mean", "flux": "flux", "u_min": "$u_{\\min}$", "v_max": "$v_{\\max}$",
           "v_min": "$v_{\\min}$", "ek": "kinetic energy", "eps": "dissipation rate",
           "umax": "$u_{\\max}$"}
    lines, last = [], None
    for k, r in frozen.items():
        if last is not None and k[0] != last:
            lines.append("\\addlinespace[1pt]")
        c = released[k]
        lines.append(f"{k[0] if k[0] != last else ''} & {obs[k[1]]} & {tex_tg(k[2])} & {r['n']} & "
                     f"{r['accepted']} & {r['distinct']} & {r['missed']} & "
                     f"{r['declined_reach'] + r['declined_unreach']} & {tex_num(r['pen_med'])} & "
                     f"{r['runs_med']:.0f} & {c['accepted']}/{c['missed']} \\\\")
        last = k[0]
    return "\n".join(lines) + "\n"


def ablation_table():
    """Table S-3.2: the released program with one check removed or changed."""
    keys = {(r["case"], r["ob"], r["target"]) for r in load("v1.2.1")}
    lines = []
    for cfg, lab in ABLATION:
        R = [r for r in load(cfg) if (r["case"], r["ob"], r["target"]) in keys]
        T = totals(R)
        where = list(dict.fromkeys(LABEL[r["case"]] for r in R if r["missed"]))
        miss = tex_int(T["miss"]) + (f" ({T['distinct_missed']})" if T["miss"] else "")
        lines.append(f"{lab} & {tex_int(T['acc'])} & {tex_int(T['distinct'])} & {miss} & "
                     f"{T['cases_miss']} & {tex_int(T['dr'])} & {tex_num(T['pen'])} & "
                     f"{tex_num(T['geo'])} & {T['runs']:.0f} & {', '.join(where) or '---'} \\\\")
    return "\n".join(lines) + "\n"


def baselines_table(base):
    """Table S-4.2: the program and the three procedures of current practice by group."""
    out = []
    for g, keys in GROUPS:
        rows = [v for k, v in base.items() if k.split("|")[0] in keys]
        for i, (key, lab) in enumerate(PROCEDURES):
            R = [v[key] for v in rows]
            out.append(f"{f'{g} ({len(rows)} cases)' if i == 0 else ''} & {lab} & "
                       f"{sum(r['reached'] for r in R)} & {sum(r['missed'] for r in R)} & "
                       f"{sum(r['declined'] for r in R) if key == 'program' else '---'} & "
                       f"{tex_num(median([r['pen'] for r in R]))} & "
                       f"{np.median([r['runs'] for r in R]):.0f} & "
                       f"{tex_num(median([r['K'] for r in R]))} \\\\")
        out.append("\\addlinespace[2pt]")
    return "\n".join(out[:-1])


def md_table():
    """Table S-4.3: LAMMPS and GROMACS double precision by size and target."""
    R = load_json("baselines_lmp_gmx.json")["cells"]
    lab = {"LMP-20": "LAMMPS 8,000", "LMP-40": "LAMMPS 64,000", "LMP-64": "LAMMPS 262,144",
           "GMX-d": "GROMACS double"}
    procs = [PROCEDURES[i] for i in (0, 1, 3, 6, 7)]
    out = ["\\begin{tabular}{ll" + "r" * len(procs) + "}", "\\toprule",
           "Code and size & Target & " + " & ".join(p[1] for p in procs) + " \\\\", "\\midrule"]
    for k, v in R.items():
        c, _, tg = k.split("|")
        cells = []
        for p, _ in procs:
            x = v[p]
            pen = "---" if not np.isfinite(x["pen"]) else f"{x['pen']:.2f}"
            cells.append(f"{x['reached']}/{x['reached'] + x['missed'] + x['declined']} ({pen})")
        out.append(f"{lab[c]} & {tex_tg(float(tg))} & " + " & ".join(cells) + " \\\\")
    return "\n".join(out + ["\\bottomrule", "\\end{tabular}"])


def transfer_table():
    """Table S-4.1: failure rate of a convergence threshold chosen once."""
    cls = {"discretisation": "d", "nested-solver": "i", "cost-only": "c"}
    par = {"n_arc": "$n_\\mathrm{arc}$", "dt": "$\\Delta t$"}
    out = ["\\begin{tabular}{lllrrrrrrrrrr}", "\\toprule",
           "& & & \\multicolumn{5}{c}{Without the adjustments} & "
           "\\multicolumn{5}{c}{With the adjustments} \\\\",
           "\\cmidrule(lr){4-8}\\cmidrule(lr){9-13}",
           "Problem & Parameter & & $f{=}0.01$ & $0.05$ & $0.1$ & $0.25$ & $0.5$ & $0.01$ & "
           "$0.05$ & $0.1$ & $0.25$ & $0.5$ \\\\",
           "\\midrule"]
    keyed = {}
    for r in load_json("transfer_floor.json"):
        keyed.setdefault((r["case"], r["param"], r["kind"]), {})[r["floor"]] = r["rates"]
    last = None
    for (case, p, kind), v in keyed.items():
        if last is not None and case != last:
            out.append("\\addlinespace[1pt]")
        out.append(f"{case if case != last else ''} & {par.get(p, p.replace('_', chr(92) + '_'))} & "
                   f"{cls[kind]} & " + " & ".join(f"{x:.0f}\\%" for x in v[False] + v[True]) + " \\\\")
        last = case
    return "\n".join(out + ["\\bottomrule", "\\end{tabular}"])


def estimate_table():
    """Table S-3.3: ratio of the estimated to the true error."""
    out = ["\\begin{tabular}{llrrrrr}", "\\toprule",
           "& & & Signed sum & \\multicolumn{3}{c}{Final estimate} \\\\", "\\cmidrule(lr){5-7}",
           "Problem & Observable & Settings & under. & 10th pct. & Median & Under. \\\\",
           "\\midrule"]
    last = None
    for r in load_json("estimate_v12.json"):
        a = r["v12"]
        out.append(f"{LABEL.get(r['case'], r['case']) if r['case'] != last else ''} & "
                   f"{OBS.get(r['ob'], r['ob'])} & {r['n']} & {r['raw']:.0f}\\% & {a[0]:.2f} & "
                   f"{a[1]:.2f} & {a[2]:.1f}\\% \\\\")
        last = r["case"]
    return "\n".join(out + ["\\bottomrule", "\\end{tabular}"])


if __name__ == "__main__":
    R = load("v1.2.1")
    base = load_json("baselines_v6.json")["cells"]
    base = {k: v for k, v in base.items() if k.split("|")[0] in LABEL}
    tables = {"tab_pooled.tex": pooled_table(R), "tab_pooled_compact.tex": compact_table(R, base),
              "tab_preregistered.tex": preregistered_table(), "tab_sup_ablation.tex": ablation_table(),
              "tab_estimate.tex": estimate_table(), "tab_transfer.tex": transfer_table(),
              "tab_baselines_full.tex": baselines_table(base), "tab_lmpbase.tex": md_table()}
    for name, text in tables.items():
        open(os.path.join(OUT, name), "w").write(text)
    T = totals(R)
    print(f"v1.2.1: {T['cases']} test cases, {T['n']} searches ({T['sf']} from unconverged starts), "
          f"{T['acc']} accepted ({T['distinct']} distinct), {T['miss']} missed, "
          f"{T['dr'] + T['du']} declined; median cost ratio {T['pen']:.2f}, runs {T['runs']:.0f}")
