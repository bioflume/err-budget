# Security policy

## Supported versions

| Version | Supported |
| ------- | --------- |
| 2.x     | Yes       |
| < 2.0   | No        |

## Reporting a vulnerability

Please report a vulnerability privately, not in a public issue. Either

- use GitHub's private reporting: on the repository page, open the
  **Security and quality** tab, click **Report a vulnerability** and fill in the form; or
- email the maintainer, Gökberk Kabacaoğlu, at <gokberk.kabacaoglu@durham.ac.uk>.

Say what the problem allows, which version you used (`errbudget version`), and how to
reproduce it.

## What to expect

You will receive an acknowledgement within a week. errbudget is research software, so the
work that follows is best effort: there is no guaranteed time for a fix, but you will be
told whether the report is accepted and kept informed of progress, and we will agree with
you when the problem is made public.

## Job files run commands by design

A job file tells errbudget what to run. The `command` field is run through the shell once
for every setting that is tried. This is how errbudget drives a solver, not a
vulnerability, and it means that a job file can do anything a script can. Treat a
job file from someone else as you would a script from them, and read it before you run
it: `errbudget check` and `errbudget run` execute it, while `errbudget check --no-run`
only validates it.

A report that a job file can run arbitrary commands is therefore out of scope. A way to
make errbudget run something that the job file does not ask for, for example by reading a
result file or a run log, is in scope.
