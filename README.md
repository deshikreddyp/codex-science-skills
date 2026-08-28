# Codex Science Skills

Private source-of-truth repository for D. Putluru's reusable Codex workflows.

## Skills

- `csf-flow-postprocess`: CSF flow rate, pressure, area deformation, and tissue displacement.
- `drug-transport-postprocess`: cross-sectional `c_net` and `c_avg` transport metrics.
- `steady-streaming-postprocess`: Eulerian, ALE, and Lagrangian steady-streaming analysis.
- `paper-style-rewrite`: scientific-paper rewriting under Gomez Lab style gates.
- `slurm-job-submission`: reliable RCAC Slurm scripts, validation, submission, and troubleshooting.
- `tex-notes`: compact `short-note-tex` briefs and explanatory `long-note-tex` scientific documentation.

## Local installation

Codex discovers global personal skills under `$HOME/.agents/skills`. Keep this repository as the canonical checkout and symlink each skill directory:

```bash
mkdir -p "$HOME/.agents/skills"
repo="$HOME/.local/share/codex-science-skills"
for skill in "$repo"/.agents/skills/*; do
  ln -s "$skill" "$HOME/.agents/skills/$(basename "$skill")"
done
```

Update the local checkout explicitly:

```bash
git -C "$HOME/.local/share/codex-science-skills" pull --ff-only
```

Codex normally detects local changes. Restart Codex if an updated skill does not appear.

## Validation

```bash
python scripts/validate_skills.py
python .agents/skills/csf-flow-postprocess/scripts/test_synthetic.py
python .agents/skills/drug-transport-postprocess/scripts/test_synthetic.py
python .agents/skills/csf-flow-postprocess/scripts/test_xdmf_end_to_end.py
python .agents/skills/drug-transport-postprocess/scripts/test_xdmf_end_to_end.py
python .agents/skills/slurm-job-submission/scripts/test_validate_slurm_script.py

# Optional on hosts with DOLFIN/FEniCS:
python .agents/skills/csf-flow-postprocess/scripts/test_fenics_h5_end_to_end.py
python .agents/skills/drug-transport-postprocess/scripts/test_fenics_h5_end_to_end.py
```

Simulation results, depot data, generated figures, credentials, plugin caches, and OpenAI-provided system skills do not belong in this repository.
