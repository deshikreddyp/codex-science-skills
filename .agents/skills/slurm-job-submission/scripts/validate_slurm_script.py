#!/usr/bin/env python3
"""Statically validate resource and reliability choices in a Slurm script."""

from __future__ import annotations

import argparse
import json
import math
import re
import shlex
import socket
from dataclasses import asdict, dataclass
from pathlib import Path


CPU_CORES = {"bell": 128, "negishi": 128, "gautschi": 192}
KEY_ALIASES = {
    "A": "account",
    "account": "account",
    "p": "partition",
    "partition": "partition",
    "N": "nodes",
    "nodes": "nodes",
    "n": "ntasks",
    "ntasks": "ntasks",
    "ntasks-per-node": "ntasks-per-node",
    "c": "cpus-per-task",
    "cpus-per-task": "cpus-per-task",
    "t": "time",
    "time": "time",
    "J": "job-name",
    "job-name": "job-name",
    "o": "output",
    "output": "output",
    "e": "error",
    "error": "error",
    "array": "array",
    "qos": "qos",
    "mem": "mem",
}
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")


@dataclass(frozen=True)
class Finding:
    level: str
    message: str
    line: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("script", type=Path)
    parser.add_argument("--cluster", choices=("auto", "bell", "negishi", "gautschi", "other"), default="auto")
    parser.add_argument("--strict", action="store_true", help="Return nonzero for warnings as well as errors")
    parser.add_argument("--json", action="store_true", help="Emit a machine-readable report")
    return parser.parse_args()


def detect_cluster(requested: str) -> str:
    if requested != "auto":
        return requested
    host = socket.getfqdn().lower()
    for name in CPU_CORES:
        if name in host:
            return name
    return "other"


def parse_directive_tokens(raw: str) -> list[tuple[str, str]]:
    try:
        tokens = shlex.split(raw, comments=False, posix=True)
    except ValueError:
        return []
    parsed: list[tuple[str, str]] = []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token.startswith("--"):
            item = token[2:]
            if "=" in item:
                key, value = item.split("=", 1)
            else:
                key = item
                if i + 1 < len(tokens) and not tokens[i + 1].startswith("-"):
                    i += 1
                    value = tokens[i]
                else:
                    value = "true"
            parsed.append((KEY_ALIASES.get(key, key), value))
        elif token.startswith("-") and len(token) == 2:
            key = token[1:]
            if i + 1 < len(tokens):
                i += 1
                parsed.append((KEY_ALIASES.get(key, key), tokens[i]))
        i += 1
    return parsed


def positive_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


def thread_value(raw: str | None, cpus_per_task: int) -> int | None:
    if raw is None:
        return None
    value = raw.strip("\"'")
    if value in ("$SLURM_CPUS_PER_TASK", "${SLURM_CPUS_PER_TASK}"):
        return cpus_per_task
    return positive_int(value)


def validate_text(text: str, script: Path, cluster: str) -> tuple[dict[str, str], list[Finding]]:
    lines = text.splitlines()
    findings: list[Finding] = []
    directives: dict[str, str] = {}
    directive_lines: dict[str, int] = {}
    first_command_seen = False

    if not lines or not lines[0].startswith("#!"):
        findings.append(Finding("error", "Missing shebang on the first line", 1))

    for lineno, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("#SBATCH"):
            if first_command_seen:
                findings.append(Finding("error", "#SBATCH directive appears after the first shell command and will be ignored", lineno))
            raw = stripped[len("#SBATCH") :].strip()
            if "$" in raw:
                findings.append(Finding("error", "Shell variables do not expand inside #SBATCH directives", lineno))
            for key, value in parse_directive_tokens(raw):
                directives[key] = value
                directive_lines[key] = lineno
        elif not stripped or stripped.startswith("#") or (lineno == 1 and stripped.startswith("#!")):
            continue
        else:
            first_command_seen = True

    placeholders = [(i, line) for i, line in enumerate(lines, 1) if re.search(r"<[A-Za-z][^>]*>", line)]
    for lineno, _ in placeholders:
        findings.append(Finding("error", "Unresolved <placeholder>", lineno))

    required = ("account", "partition", "time", "job-name", "output", "error")
    for key in required:
        if key not in directives:
            findings.append(Finding("error", f"Missing #SBATCH --{key}"))

    if "nodes" not in directives:
        findings.append(Finding("warning", "Specify #SBATCH --nodes explicitly"))
    if "ntasks" not in directives and "ntasks-per-node" not in directives:
        findings.append(Finding("error", "Specify --ntasks or --ntasks-per-node"))
    if "cpus-per-task" not in directives:
        findings.append(Finding("warning", "Specify --cpus-per-task explicitly"))

    nodes = positive_int(directives.get("nodes")) or 1
    ntasks = positive_int(directives.get("ntasks"))
    tasks_per_node = positive_int(directives.get("ntasks-per-node"))
    cpus_per_task = positive_int(directives.get("cpus-per-task")) or 1

    for key in ("nodes", "ntasks", "ntasks-per-node", "cpus-per-task"):
        if key in directives and positive_int(directives[key]) is None:
            findings.append(Finding("error", f"--{key} must be a positive integer", directive_lines.get(key)))

    if ntasks is not None and tasks_per_node is not None and ntasks != nodes * tasks_per_node:
        findings.append(Finding("error", f"--ntasks={ntasks} conflicts with --nodes={nodes} and --ntasks-per-node={tasks_per_node}"))
    if tasks_per_node is None and ntasks is not None:
        tasks_per_node = math.ceil(ntasks / nodes)
        if nodes > 1:
            findings.append(Finding("warning", "Multi-node placement is clearer with --ntasks-per-node"))

    cores = CPU_CORES.get(cluster)
    if cores is not None and tasks_per_node is not None:
        requested_per_node = tasks_per_node * cpus_per_task
        if requested_per_node > cores:
            findings.append(Finding("error", f"Requests {requested_per_node} CPUs/node but standard {cluster} CPU nodes have {cores}"))
        elif requested_per_node == cores:
            findings.append(Finding("info", f"Requests one full standard {cluster} CPU node ({cores} CPUs/node)"))

    body = "\n".join(lines)
    if ("source activate" in body or "conda activate" in body or "pipefail" in body) and lines:
        if "bash" not in lines[0]:
            findings.append(Finding("error", "Conda activation or pipefail requires a Bash shebang", 1))
    if not re.search(r"^\s*set\s+-[^\n]*e", body, re.MULTILINE):
        findings.append(Finding("warning", "Use `set -e` (normally `set -euo pipefail`) for fail-fast behavior"))
    if "module --force purge" not in body and "module purge" not in body:
        findings.append(Finding("warning", "Purge inherited modules before loading the job environment"))
    if re.search(r"\bcd\s+\$SLURM_SUBMIT_DIR\b", body):
        findings.append(Finding("warning", "Quote the submission directory: cd -- \"$SLURM_SUBMIT_DIR\""))

    exports = {
        match.group(1): match.group(2).strip()
        for match in re.finditer(r"^\s*export\s+([A-Za-z_][A-Za-z0-9_]*)=(.+?)\s*$", body, re.MULTILINE)
    }
    rank_count = tasks_per_node or ntasks or 1
    if rank_count > 1 and cpus_per_task == 1:
        for variable in THREAD_VARS:
            if variable not in exports:
                findings.append(Finding("warning", f"Pure MPI job does not cap {variable}=1"))
            else:
                parsed_threads = thread_value(exports[variable], cpus_per_task)
                if parsed_threads is None:
                    findings.append(Finding("warning", f"Cannot verify that pure MPI setting {variable} equals 1"))
                elif parsed_threads != 1:
                    findings.append(Finding("error", f"Pure MPI job sets {variable} above 1 and may oversubscribe CPUs"))
        launch_lines = [line for line in lines if re.search(r"(^|\s)(srun|mpirun|mpiexec)(\s|$)", line) and not line.lstrip().startswith("#")]
        if launch_lines and "--cpu-bind" not in body and "--bind-to" not in body:
            findings.append(Finding("warning", "Pure MPI launcher does not explicitly bind ranks to cores"))
        if any(re.search(r"(^|\s)srun(\s|$)", line) for line in launch_lines) and "--kill-on-bad-exit" not in body:
            findings.append(Finding("warning", "Pure MPI srun launcher does not request --kill-on-bad-exit=1"))
    elif cpus_per_task > 1:
        omp = thread_value(exports.get("OMP_NUM_THREADS"), cpus_per_task)
        if "OMP_NUM_THREADS" not in exports:
            findings.append(Finding("warning", "Multi-CPU task does not set OMP_NUM_THREADS; choose 1 for worker pools or SLURM_CPUS_PER_TASK for threaded code"))
        elif omp not in (1, cpus_per_task):
            findings.append(Finding("warning", f"OMP_NUM_THREADS does not match either 1 or --cpus-per-task={cpus_per_task}"))
        for variable in THREAD_VARS[1:]:
            if variable not in exports:
                findings.append(Finding("warning", f"Multi-CPU Python job does not bound {variable}; set it deliberately to prevent nested threading"))

    for key in ("output", "error"):
        value = directives.get(key)
        if value and not value.startswith("/"):
            findings.append(Finding("warning", f"--{key} is relative; confirm the submission directory is intentional", directive_lines.get(key)))
        if value and value.startswith("/"):
            static_parent = Path(re.sub(r"%[A-Za-z]", "token", value)).parent
            if not static_parent.exists():
                findings.append(Finding("warning", f"Parent directory for --{key} does not currently exist: {static_parent}", directive_lines.get(key)))

    for lineno, line in enumerate(lines, 1):
        if re.search(r"\b(?:srun|mpirun|mpiexec)\b.*(?:^|\s)>+\s*[^&]", line):
            findings.append(Finding("warning", "Launcher stdout is redirected even though Slurm output is configured", lineno))

    return directives, findings


def main() -> None:
    args = parse_args()
    if not args.script.is_file():
        raise SystemExit(f"Script not found: {args.script}")
    cluster = detect_cluster(args.cluster)
    directives, findings = validate_text(args.script.read_text(encoding="utf-8"), args.script, cluster)
    errors = sum(item.level == "error" for item in findings)
    warnings = sum(item.level == "warning" for item in findings)
    report = {
        "script": str(args.script.resolve()),
        "cluster": cluster,
        "directives": directives,
        "findings": [asdict(item) for item in findings],
        "summary": {"errors": errors, "warnings": warnings},
    }
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        for item in findings:
            where = f" line {item.line}" if item.line is not None else ""
            print(f"{item.level.upper()}{where}: {item.message}")
        if not findings:
            print("OK: no findings")
        print(f"Summary: {errors} error(s), {warnings} warning(s); cluster={cluster}")
    raise SystemExit(1 if errors or (args.strict and warnings) else 0)


if __name__ == "__main__":
    main()
