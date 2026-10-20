# Example 4: a solver that is a command

`solver.py` stands in for a simulation code that takes its numerical parameters on the
command line and prints a result. `job.yaml` tells `errbudget` how to run it, which
parameters to tune and what accuracy to ask for.

```bash
pip install "errbudget[yaml]"
errbudget check job.yaml   # validates the file and runs the starting setting once
errbudget run job.yaml     # tunes (about 60 runs, a few seconds), prints the report
errbudget runs job.yaml    # the run log
```

The exact solution is `cos(1) + 0`, so you can check the advice: the recommended setting
is about 15x cheaper than the cautious start and still within 1e-4.
