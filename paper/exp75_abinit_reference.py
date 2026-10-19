"""Revision 4: the uncertainty of the ABINIT reference that the tensor is graded against.

Until revision 4, cases.py used the uncertainty measured in revision 2 for a 30-Ha
reference (figures/sec83_reference.txt: 1.22e-4 in the energy, 4.45e-3 in the force),
although the tensor (abinit_tensor_v11.pkl) is graded against the 48-Ha reference
(abinit_silicon.py, reference(): ecut 48 Ha, ngkpt 8, toldfe 1e-14, nline 8).  This
measures the uncertainty of the 48-Ha reference with the recipe used for Quantum
ESPRESSO (exp53_qe.py, hamilton/J3_qe): the change from one cutoff step below the
reference plus the change from 8 to 10 k-points at a lower cutoff,

    u = (|X(40 Ha, 8) - X(48 Ha, 8)| + |X(36 Ha, 10) - X(36 Ha, 8)|) / |X(48 Ha, 8)|,

and, for comparison, the rule of exp45_abinit_reference.py (the larger of the changes
42 -> 48 Ha and 8 -> 10 k-points at 48 Ha).  It checks that the 48-Ha run reproduces the
stored reference and writes data/abinit_reference_v6.json, which cases.py then uses
(set KNOBS_ABINIT_UNC=rev2 to keep the old value).

    python exp75_abinit_reference.py        # 6 ABINIT runs, about an hour on 2 cores
"""
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from abinit_silicon import AbinitSilicon

RUNS = [("ref", {}), ("e40", dict(ecut=40.0)), ("e42", dict(ecut=42.0)),
        ("k8", dict(ecut=36.0)), ("k10", dict(ecut=36.0, ngkpt=10)), ("ref_k10", dict(ngkpt=10))]
OBS = ("energy", "force")


def main():
    p = AbinitSilicon()
    base = p.reference()
    _, STORED = pickle.load(open(os.path.join(HERE, "data", "abinit_tensor_v11.pkl"), "rb"))
    vals, runs = {}, {}
    for lab, upd in RUNS:
        th = dict(base, **upd)
        obs, cost, wall, info = p.run(th)
        if not all(k in obs and obs[k] == obs[k] for k in OBS):
            raise SystemExit(f"ABINIT run {lab} {th} failed: {info}")
        vals[lab] = {k: float(obs[k]) for k in OBS}
        runs[lab] = dict(theta=th, obs=vals[lab], cpu_s=cost, wall_s=wall,
                         converged=bool(obs.get("converged", True)))
        print(f"{lab:<8} {th}: E = {obs['energy']:.10f} Ha, |F| = {obs['force']:.8e} Ha/Bohr "
              f"(cpu {cost:.0f} s, converged {obs.get('converged')})", flush=True)
    R = vals["ref"]
    same = all(abs(R[o] - STORED[o]) <= 1e-8 * abs(STORED[o]) for o in OBS)
    print(("REFERENCE REPRODUCED" if same else "REFERENCE NOT REPRODUCED") +
          f": stored {dict((o, STORED[o]) for o in OBS)}, this run {R}")
    unc, unc_max_rule, parts = {}, {}, {}
    for o in OBS:
        de = abs(vals["e40"][o] - R[o]) / abs(R[o])
        dk = abs(vals["k10"][o] - vals["k8"][o]) / abs(R[o])
        unc[o] = de + dk
        de42 = abs(vals["e42"][o] - R[o]) / abs(R[o])
        dk48 = abs(vals["ref_k10"][o] - R[o]) / abs(R[o])
        unc_max_rule[o] = max(de42, dk48)
        parts[o] = dict(cutoff_40_48=de, kpoints_8_10_at_36=dk, cutoff_42_48=de42,
                        kpoints_8_10_at_48=dk48)
        print(f"{o:<6}: 40 -> 48 Ha {de:.2e}, 8 -> 10 k-points at 36 Ha {dk:.2e}; "
              f"uncertainty (QE recipe) {unc[o]:.2e};  exp45 rule (42 -> 48 Ha {de42:.2e}, "
              f"8 -> 10 k-points at 48 Ha {dk48:.2e}) {unc_max_rule[o]:.2e};  "
              f"revision 2 value {dict(energy=1.22e-4, force=4.45e-3)[o]:.2e}")
    for o, tgs in (("energy", [1e-2, 3e-3, 1e-3]), ("force", [1e-1, 3e-2, 1e-2])):
        for tg in tgs:
            print(f"  {o:<6} target {tg:.0e}: {tg / unc[o]:9.1f} times the uncertainty -> "
                  f"{'admitted' if tg >= 10 * unc[o] else 'not admitted'} "
                  f"(revision 2: {'admitted' if tg >= 10 * dict(energy=1.22e-4, force=4.45e-3)[o] else 'not admitted'})")
    out = dict(ref=R, stored_ref={o: float(STORED[o]) for o in OBS}, reproduced=same, unc=unc,
               unc_exp45_rule=unc_max_rule, parts=parts, runs=runs,
               recipe="(|X(40,8)-X(48,8)| + |X(36,10)-X(36,8)|)/|X(48,8)|")
    path = os.path.join(HERE, "data", "abinit_reference_v6.json")
    json.dump(out, open(path, "w"), indent=1)
    print("written", path)


if __name__ == "__main__":
    main()
