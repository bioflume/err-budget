# Data kept outside the repository

The repository holds the code, the stored solver runs that the replays read (`paper/data/`, `hamilton/results/`) and
the small summaries of the cluster studies (`studies/*/summary/`). The raw outputs below are too large for git. Each
is listed with its SHA-256, so that an archived copy can be checked. They are to be deposited in one Zenodo record
(DOI to be added).

Laptop paths are on the author's Mac. Hamilton paths are under `/nobackup/grbj73/`.

## Hamilton jobs J1–J5 (`hamilton/`)

| File | Size | SHA-256 | Laptop | Contents |
|---|---:|---|---|---|
| `hamilton_results.tar.gz` | 46.8 MB | `7a98b2b6247a2199cba2c9abbdb2b543a693b9b484c7378c85d3b0dc930858cd` | `~/Desktop/knobs-runs/` | J1–J5, 31 October, with `results/raw.tar.gz` |
| `hamilton_results.tar.gz` | 46.9 MB | `dd029a8bb54ec55e9a1595c68f45cd9bd45cc2f4e1e0b7718f9effc1b9db9ad2` | `~/Desktop/knobs-runs/rerun/results/` | J2 and J4 reruns, 3 November, first pack |
| `hamilton_results.tar.gz` | 47.5 MB | `b142724a087c292e3bbc90eff42fce975d0858b060519bf2f95295916c15bdc0` | `~/Desktop/knobs-runs/rerun/results2/` | J2 and J4 reruns, 3 November, final pack |
| `errbudget_v6_results.tar.gz` | 0.4 MB | `d1bfea0235899eb599f492561dafdb9081edf6a6528b8dd41cdbbafe88c48746` | `~/Desktop/knobs-runs/analysis_v3/` | `hamilton/results` and the replay outputs as collected on 3 November; the tracked files of `hamilton/results` come from it |

On Hamilton the jobs ran in `/nobackup/grbj73/errbudget/` (`hamilton/config.sh`).

## Reynolds-number reuse study (`studies/reuse2d1/`)

| File | Size | SHA-256 | Contents |
|---|---:|---|---|
| `dist/reuse2d1-batch.tar.gz` | 27.5 kB | `95aa1939e9eb5a37bd6590d1d48eb1907acb4ae991892edc504441baa4dc4163` | the code sent to Hamilton (project commit 530c4b6, code 9b7482e) |
| `results/reuse2d1-20261115-1111.tar.gz` | 2.5 MB | `b13f990859632408ddee8b6502ea527f5f98fbe8e074d1bc14cbe7ec2a4380fd` | every run of the grid and the references (43 MB unpacked) |

Laptop: `~/Desktop/knobs-runs/finalResults/reuse2d1/{dist,results}/`. Hamilton: `/nobackup/grbj73/reuse2d1/`.
`summary/reuse2d1-20261115-1111/` holds the analysis output of this pack (`res.txt`, `reuse_summary.json`).

## Airfoil study (`studies/airfoil/`)

| File | Size | SHA-256 | Contents |
|---|---:|---|---|
| `dist/airfoil-batch.tar.gz` | 141 kB | `360b1b94248c067af39f35d2fe95f17a1fdc1916773fa53acb1516a070dbfce0` | the last code sent to Hamilton (project commit a8d5701, amendment A5) |
| `results/airfoil-check-A2.tar.gz` | 115 kB | `a81d28bd2a7254a3cf22a22bdc506bd85b5ef8d88a26ca6a827cb800ebab9912` | check job after amendment A2 |
| `results/airfoil-diag-A2.tar.gz` | 58 kB | `a328d61fc68dccf6dadc9134dc8973c5ded0bb4ec8136a379acf71626c0f3a53` | diagnostic runs after A2 |
| `results/check-runs.tar.gz` | 55 kB | `93c57fb5c94ad0630405de822209471a8c5dcee8882e4ef7e21ee4200e9ee9ff` | the two check runs of the convergence map |
| `results/airfoil-convmap-20261111-0910.tar.gz` | 34.0 MB | `81718623256aa46a00055568941c790c47da931bc03282ad2047421ef6d5a1fc` | convergence map (A3) |
| `results/airfoil-confirm-20261112-0927.tar.gz` | 36.8 MB | `c61f8a1401dc8ac92f1dcadf65ed99710df854fd7bdb9cba09e5ab6505c12ac9` | confirmation map (A4) |
| `results/airfoil-settle-20261113-1458.tar.gz` | 37.9 MB | `f5d2112a18479fb6138198add802d34ec40949b096d364a88225d4ab6c4adfd0` | settling measurement |
| `results/airfoil-stage0-20261113-2123.tar.gz` | 38.6 MB | `1b5685e69bdcfceab7228e2dcd5501aedca1c930f744e0d208c8bebe69cbe898` | Stage 0, 25 of 26 runs |
| `results/airfoil-stage0-20261113-2137.tar.gz` | 38.6 MB | `fa59fb6ddaf09c124077319747a63ed090dfe76dd3f0ca94132c0085c4698253` | Stage 0, all 26 runs |
| `results/airfoil-stage1-20261114-1548.tar.gz` | 52.5 MB | `026d5f00ba872ab78ed0a6ab7e800dd52440202f79b06bf45602d8566a84477a` | Stage 1, 546 runs |

Laptop: `~/Desktop/knobs-runs/finalResults/airfoil/{dist,results}/`; `results/local/` there holds the local Docker
runs of the 10 November diagnosis. Hamilton: `/nobackup/grbj73/airfoil/` (packs in `packs/`). Unpacked, the results
take 6.4 GB. `summary/` holds the analysis output of the Stage 0, Stage 1, convergence-map and confirmation-map
packs, each in a folder named after its pack.

## Not archived

OpenFOAM case directories are not kept: each run's driver stores its observables, iteration counts and logs in the
packs above. The PRE comparison downloads the authors' code (`studies/pre_comparison/fetch_pre.sh`) and keeps
nothing of it.
