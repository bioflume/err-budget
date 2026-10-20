"""Post-hoc variant examined after the held-out replay of H2 (not part of the frozen
program).  In version 1.2 a parameter takes part in the joint-refinement check only if
its one-step change exceeds three times the resolution (0.3 eps), so on H2 neither
parameter qualified and the check was skipped at the setting that missed.  The variant
lowers this to the resolution itself (0.1 eps).  Import this module to patch
errbudget._core in the current process.
"""
import inspect
import textwrap

import errbudget._core as K

src = textwrap.dedent(inspect.getsource(K.LocalBudget.joint_correct))
assert "abs(self.e[k][j]) <= 3.0 * self.resolution" in src
src = src.replace("abs(self.e[k][j]) <= 3.0 * self.resolution", "abs(self.e[k][j]) <= self.resolution")
ns = {}
exec(compile(src, "<posthoc_v121>", "exec"), K.__dict__, ns)
K.LocalBudget.joint_correct = ns["joint_correct"]
