#!/usr/bin/env python3
"""Self-tests for validate_slurm_script.py."""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

from validate_slurm_script import validate_text


GOOD_MPI = """#!/bin/bash -l
#SBATCH --job-name=flow
#SBATCH --account=test-account
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=/tmp/%x.o%j
#SBATCH --error=/tmp/%x.e%j
set -euo pipefail
module --force purge
module load anaconda
source activate fenicsproject
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
cd -- "$SLURM_SUBMIT_DIR"
srun --mpi=pmix_v3 --ntasks="$SLURM_NTASKS" --cpu-bind=cores --kill-on-bad-exit=1 python3 -u solve.py
"""


def levels(text: str) -> list[tuple[str, str]]:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "job.sbatch"
        path.write_text(text, encoding="utf-8")
        _, findings = validate_text(text, path, "gautschi")
    return [(item.level, item.message) for item in findings]


def main() -> None:
    good = levels(GOOD_MPI)
    assert not [item for item in good if item[0] == "error"], good

    oversubscribed = GOOD_MPI.replace("--cpus-per-task=1", "--cpus-per-task=2")
    over = levels(oversubscribed)
    assert any(level == "error" and "384 CPUs/node" in message for level, message in over), over

    uncapped = GOOD_MPI.replace("export OPENBLAS_NUM_THREADS=1\n", "")
    missing = levels(uncapped)
    assert any(level == "warning" and "OPENBLAS_NUM_THREADS" in message for level, message in missing), missing

    late = GOOD_MPI.replace("#SBATCH --error=/tmp/%x.e%j\n", "") + "#SBATCH --error=/tmp/%x.e%j\n"
    late_findings = levels(late)
    assert any(level == "error" and "will be ignored" in message for level, message in late_findings), late_findings

    variable = GOOD_MPI.replace("--nodes=1", "--nodes=$NODES")
    variable_findings = levels(variable)
    assert any(level == "error" and "do not expand" in message for level, message in variable_findings), variable_findings

    assets = Path(__file__).resolve().parents[1] / "assets"
    with tempfile.TemporaryDirectory() as tmp:
        for source in sorted(assets.glob("*.sbatch")):
            rendered = re.sub(r"<[A-Za-z][^>]*>", "placeholder", source.read_text(encoding="utf-8"))
            target = Path(tmp) / source.name
            target.write_text(rendered, encoding="utf-8")
            subprocess.run(["bash", "-n", str(target)], check=True)

    print("Slurm validator tests passed")


if __name__ == "__main__":
    main()
