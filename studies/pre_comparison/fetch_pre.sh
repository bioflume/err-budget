#!/bin/sh
# The authors' PRE code at the commit used here, and the Python environment it needs.
set -e
cd "$(dirname "$0")"
mkdir -p ext
[ -d ext/sparse-probabilistic-richardson-extrapolation ] || \
  git clone https://github.com/NewcastleRSE/sparse-probabilistic-richardson-extrapolation.git ext/sparse-probabilistic-richardson-extrapolation
git -C ext/sparse-probabilistic-richardson-extrapolation checkout 3918107
python3 -m venv .venv
.venv/bin/pip install numpy==2.3.1 scipy==1.18.1 jax==0.7.2 tqdm==4.67.1 matplotlib==3.10.3 \
  sympy==1.14.0 pandas==2.3.2 scikit-learn==1.7.0
