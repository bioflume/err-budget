# Use of AI in this repository

**Responsibility.** Gökberk Kabacaoğlu (Durham University) is responsible for this software
and for every result it produces. An AI coding agent wrote most of the code under the author's
direction; the agent is not an author and carries no responsibility.

## Who wrote what

- **Model.** Claude Opus 5.5 (Anthropic), run as Claude Code. Commits it made or co-made
  end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Count on `main`, 26 November 2026.** 20 of the 30 commits carry that line. The other
  ten are README edits and merges by the author, the initial commit, and three bulk imports
  of repositories prepared elsewhere: 27f6221 (knobs 1.3.0), 015bf6b (job-file fixes from
  review) and 2fc31e5 (the refactored package). The author wrote these commits (the
  author's statement, 27 November 2026). The agent did the audit of 26 November and engine
  1.2.2 in ten commits that carry the line; the author copied their result into this
  repository as one commit, 3b57805 ("Release version updated", 27 November 2026), which
  does not carry it. Every later commit made by the agent does.
- **The engine.** `src/errbudget/_core.py` is engine 1.2.2. Its configuration `v1.2.1` runs
  the engine that produced the article's results on the branch `paper`, search for search
  (`TRUST.md`). The corrections of 1.2.2 were made at the author's request of
  27 November 2026, after an independent review.

## What the author owns

The author, not the agent, owns:

- the requirements: the method and its rules, as stated in the article;
- the acceptance tests: `tests/test_spec_v121.py`, `tests/test_claims.py` and
  `tests/test_v122.py` (written by the agent, accepted by the author on 27 November 2026);
- every decision to accept a change.

These files are listed in `PROTECTED`. CI checks them against `PROTECTED.sha256`, and
`AGENTS.md` forbids agents to change them.

## Evidence that does not rest on trusting the author

`TRUST.md` lists the checks a reader can run and their results.

## Session logs

The agent's session logs (prompts, tool calls, its own reports) are not in this repository.
The 20 commits made in agent sessions name their session in a `Claude-Session:` line, three
sessions in all. Those logs were not kept, so the record of how the code was written is the
git history alone.
