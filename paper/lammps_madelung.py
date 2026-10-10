"""A real third-party solver: LAMMPS electrostatics, with analytic ground truth.

The rocksalt lattice is a simple cubic arrangement of alternating unit charges.
Its electrostatic energy is known in closed form through the Madelung constant,

    E = - (N/2) * alpha * k_e * e^2 / r0 ,   alpha = 1.7475645946...

so the error in LAMMPS's computed Coulomb energy is exactly knowable.  The
parameters exercised are LAMMPS's own, documented ones:

    pair_style lj/cut/coul/long <cutoff>   real-space cutoff   discretisation
    kspace_style pppm <accuracy>           reciprocal accuracy discretisation
    kspace_modify order <n>                interpolation order discretisation
    neigh_modify every <n> delay <n>       neighbour rebuild   cost-only

LAMMPS documents `accuracy` as "the relative RMS error in per-atom forces", used
"in conjunction with the pairwise cutoff to determine ... the grid size for style
pppm".  That is a per-code, hand-derived instance of exactly the coupling this
paper measures automatically, which makes it a fair basis for comparison.
"""
from __future__ import annotations

import time

import numpy as np

MADELUNG = 1.7475645946331821906362120355443974
KE = 14.399645478425668          # e^2 / (4 pi eps0) in eV Angstrom (metal units)


class LammpsMadelung:
    """Single-point electrostatic energy of a rocksalt lattice, via LAMMPS."""

    def __init__(self, n_cells: int = 20, r0: float = 2.82, n_eval: int = 5):
        self.n, self.r0 = n_cells, r0          # r0 = NaCl nearest-neighbour, Angstrom
        self.natoms = n_cells ** 3
        self.n_eval = n_eval                   # force evaluations timed per run

    def exact_observables(self) -> dict:
        e = -(self.natoms / 2.0) * MADELUNG * KE / self.r0
        return {"coulomb_energy": e, "energy_per_ion": e / self.natoms}

    def run(self, theta: dict):
        import lammps

        cutoff = float(theta["cutoff"])
        acc = float(theta["kspace_accuracy"])
        order = int(theta["kspace_order"])
        every = int(theta.get("neigh_every", 1))

        t0 = time.perf_counter()
        L = lammps.lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"])
        c = L.command
        c("units metal")
        c("atom_style charge")
        c("boundary p p p")
        c(f"lattice sc {self.r0}")
        c(f"region box block 0 {self.n} 0 {self.n} 0 {self.n}")
        c("create_box 2 box")
        c("create_atoms 1 box")
        c("mass 1 22.99")
        c("mass 2 35.45")

        # alternate the charge with the parity of the lattice site
        x = L.numpy.extract_atom("x")
        idx = np.rint(x / self.r0).astype(int)
        parity = (idx.sum(axis=1) % 2).astype(bool)
        q = L.numpy.extract_atom("q")
        q[:] = np.where(parity, -1.0, 1.0)
        t = L.numpy.extract_atom("type")
        t[:] = np.where(parity, 2, 1)

        # A combination of small cutoff, tight accuracy and low interpolation
        # order can demand a reciprocal-space grid too large to allocate.  That
        # is a property of the parameter space, not an error: the configuration
        # is unavailable on this machine, and is reported as such.
        c(f"pair_style lj/cut/coul/long 1.0 {cutoff}")
        c("pair_coeff * * 0.0 1.0")            # electrostatics only
        c(f"kspace_style pppm {acc}")
        c(f"kspace_modify order {order}")
        c(f"neigh_modify every {every} delay 0 check no")
        try:
            c("run 0 post no")                 # setup and first evaluation
        except Exception as exc:
            L.close()
            return ({"coulomb_energy": float("nan"),
                     "energy_per_ion": float("nan")},
                    float("inf"), time.perf_counter() - t0,
                    {"unavailable": str(exc)[:120]})
        setup = time.perf_counter() - t0

        # Cost is the per-evaluation compute, not the one-off setup: in
        # production the setup is amortised over many steps.
        t1 = time.perf_counter()
        for _ in range(self.n_eval):
            c("run 0 post no")
        per_eval = (time.perf_counter() - t1) / self.n_eval

        ecoul = L.get_thermo("ecoul")
        elong = L.get_thermo("elong")
        L.close()

        e_tot = float(ecoul + elong)
        return ({"coulomb_energy": e_tot, "energy_per_ion": e_tot / self.natoms},
                float(per_eval), setup, {"setup_s": setup})

    def runner(self, observable: str):
        def f(theta):
            obs, cost, _, _ = self.run(theta)
            return obs[observable], cost
        return f


if __name__ == "__main__":
    p = LammpsMadelung(n_cells=10)
    ex = p.exact_observables()
    print(f"rocksalt {p.n}^3 = {p.natoms} ions, r0 = {p.r0} A")
    print(f"Madelung energy      : {ex['coulomb_energy']:.6f} eV "
          f"({ex['energy_per_ion']:.6f} eV/ion)")
    print()
    print("cutoff  accuracy  order   LAMMPS energy      rel. error    wall")
    for cut in [6.0, 9.0, 12.0]:
        for acc in [1e-2, 1e-4, 1e-6]:
            obs, cost, w, _ = p.run(dict(cutoff=cut, kspace_accuracy=acc,
                                         kspace_order=5))
            rel = abs(obs["coulomb_energy"] - ex["coulomb_energy"]) / abs(
                ex["coulomb_energy"])
            print(f"{cut:5.1f}   {acc:.0e}     5    {obs['coulomb_energy']:14.6f}"
                  f"    {rel:.3e}    {w:.2f}s")
