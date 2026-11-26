# Use of AI in this repository

**Responsibility.** Gökberk Kabacaoğlu (Durham University) is responsible for this software
and for every result it produces. An AI coding agent wrote most of the code and the
analyses under the author's direction; the agent is not an author and carries no
responsibility.

## Who wrote what

- **Model.** Claude Opus 5.5 (Anthropic), run as Claude Code. Commits it made or co-made
  end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; one, a969925
  (revision 4 analyses), has Claude as its author.
- **Count on `paper`, 26 November 2026.** 39 of the 48 commits carry that line. The other
  nine are README edits by the author (7e908ec, 94a719a, cf64223), the initial commit, and
  imports or restructurings of material prepared elsewhere: db92080 (revision 4 results),
  f2ee6e5 (restructure for revision 4), 774d447 (merge of the knobs-paper repository) and
  150570a (the refactored branch). The trailer does not show who wrote the imported
  material. [Author to confirm: whether the agent wrote it.]
- **Studies imported as run.** `studies/` was imported on 15 November 2026 from a project
  repository that is not available. The commits that `paper/README.md` cites as the
  registrations of the three studies (airfoil da6e9e0 and amendments A1-A6, reuse2d1
  8cb358c, PRE e1d1bce) are not in this repository or in any other repository of the
  author, so their dates cannot be checked here. `TRUST.md` says what evidence remains.

## What the author owns

The author, not the agent, owns:

- the requirements: the method and its rules, as stated in the article, and the
  pre-registration (`PREREGISTRATION.md`, the study protocols);
- the acceptance tests: `tests/test_spec_v121.py` and `tests/test_claims.py`;
- the record of the results: the stored runs and printed outputs;
- every decision to accept a change or a result.

These files are listed in `PROTECTED`. CI checks them against `PROTECTED.sha256`, and
`AGENTS.md` forbids agents to change them. [Author to confirm: that the author, not an
agent, submitted the Hamilton jobs, as the runbook `hamilton/README.md` is written for.]

## Evidence that does not rest on trusting the author

`TRUST.md` lists the checks a reader can run and their results.

## Not recorded

The agent's session logs (prompts, tool calls, its own reports) were not kept, so the
record of how the code was written is the git history alone.
