"""Write ro-crate-metadata.json: every stored result, its SHA-256 and the script that writes it.

The crate (RO-Crate 1.1) describes the stored runs and printed outputs of the article
(paper/data/, paper/figures/, hamilton/results/, studies/*/summary/). Each file gets its
size and SHA-256. Where a row of the tables in paper/README.md lists a file under "Writes",
the crate records a CreateAction from each script of that row to the file. A row can list
several scripts, so a link says "made by one of these scripts", not which one.

    python tools/make_ro_crate.py      # rewrites ro-crate-metadata.json
"""
import glob
import hashlib
import json
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATTERNS = ["paper/data/**", "paper/figures/**", "hamilton/results/**", "studies/*/summary/**",
            "studies/pre_comparison/results/**"]


def tracked():
    r = subprocess.run(["git", "-C", ROOT, "ls-files", "-z"], capture_output=True, check=True)
    return set(r.stdout.decode().split("\0")) - {""}


def writers():
    """{file: [scripts]} from the tables of paper/README.md (columns Script and Writes)."""
    out = {}
    for line in open(os.path.join(ROOT, "paper", "README.md"), encoding="utf-8"):
        cols = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cols) != 4 or cols[0] in ("Result", "---"):
            continue
        scripts = re.findall(r"`([^`]+\.py)`", cols[1])
        for f in re.findall(r"`([^`]+)`", cols[3]):
            f = f if "/" in f.lstrip("./") else f"figures/{f}" if f.endswith((".txt", ".tex", ".pdf")) else f
            p = os.path.normpath(os.path.join("paper", f))
            out.setdefault(p, [])
            for s in scripts:
                s = os.path.normpath(os.path.join("paper", s))
                if s not in out[p]:
                    out[p].append(s)
    return out


def main():
    files = tracked()
    commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True,
                            text=True, check=True).stdout.strip()
    data = sorted(f for p in PATTERNS for f in glob.glob(os.path.join(ROOT, p), recursive=True)
                  if os.path.isfile(f) and os.path.relpath(f, ROOT) in files)
    data = [os.path.relpath(f, ROOT) for f in data]
    made_by = writers()
    graph = [
        {"@id": "ro-crate-metadata.json", "@type": "CreativeWork",
         "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"}, "about": {"@id": "./"}},
        {"@id": "./", "@type": "Dataset",
         "name": "errbudget: the stored runs and results of the article",
         "description": "Stored solver runs, replay outputs and cluster summaries that reproduce "
                        "the article and its supplement; see paper/README.md and DATA.md.",
         "license": {"@id": "https://spdx.org/licenses/MIT"},
         "author": {"@id": "#author"},
         "version": f"git {commit[:7]}",
         "hasPart": [{"@id": f} for f in data]},
        {"@id": "#author", "@type": "Person", "name": "Gökberk Kabacaoğlu",
         "affiliation": "Department of Computer Science, Durham University"},
    ]
    scripts = set()
    for f in data:
        p = os.path.join(ROOT, f)
        graph.append({"@id": f, "@type": "File", "contentSize": os.path.getsize(p),
                      "sha256": hashlib.sha256(open(p, "rb").read()).hexdigest()})
        for s in made_by.get(f, []):
            if s in files:
                scripts.add(s)
                graph.append({"@id": f"#write-{f}-{s}", "@type": "CreateAction",
                              "instrument": {"@id": s}, "result": {"@id": f},
                              "description": f"listed with {s} in one row of paper/README.md "
                                             f"(commit {commit[:7]})"})
    for s in sorted(scripts):
        graph.append({"@id": s, "@type": ["File", "SoftwareSourceCode"],
                      "programmingLanguage": "Python"})
    json.dump({"@context": "https://w3id.org/ro/crate/1.1/context", "@graph": graph},
              open(os.path.join(ROOT, "ro-crate-metadata.json"), "w"), indent=1,
              ensure_ascii=False)
    linked = sum(1 for f in data if any(s in files for s in made_by.get(f, [])))
    print(f"{len(data)} files, {linked} linked to the scripts of their row in paper/README.md")


if __name__ == "__main__":
    main()
