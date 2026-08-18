#!/usr/bin/env python3
"""Validate every skill's required files, frontmatter, and Python syntax."""

from __future__ import annotations

import compileall
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ROOT / ".agents" / "skills"
NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")


def validate(skill: Path) -> list[str]:
    errors: list[str] = []
    document = skill / "SKILL.md"
    if not document.is_file():
        return [f"{skill.name}: missing SKILL.md"]
    text = document.read_text(encoding="utf-8")
    match = re.match(r"^---\n(.*?)\n---\n", text, re.DOTALL)
    if not match:
        return [f"{skill.name}: invalid YAML frontmatter delimiters"]
    try:
        metadata = yaml.safe_load(match.group(1))
    except yaml.YAMLError as exc:
        return [f"{skill.name}: invalid YAML: {exc}"]
    if set(metadata or {}) != {"name", "description"}:
        errors.append(f"{skill.name}: frontmatter must contain only name and description")
    name = (metadata or {}).get("name")
    if name != skill.name or not isinstance(name, str) or not NAME_RE.fullmatch(name):
        errors.append(f"{skill.name}: invalid or mismatched name {name!r}")
    description = (metadata or {}).get("description")
    if not isinstance(description, str) or len(description.strip()) < 30:
        errors.append(f"{skill.name}: description is missing or too short")
    if not (skill / "agents" / "openai.yaml").is_file():
        errors.append(f"{skill.name}: missing recommended agents/openai.yaml")
    if (skill / "scripts").is_dir() and not compileall.compile_dir(skill / "scripts", quiet=1):
        errors.append(f"{skill.name}: Python compilation failed")
    return errors


def main() -> int:
    if not SKILLS.is_dir():
        print(f"Missing skills directory: {SKILLS}", file=sys.stderr)
        return 1
    errors = [error for skill in sorted(SKILLS.iterdir()) if skill.is_dir() for error in validate(skill)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Validated {sum(path.is_dir() for path in SKILLS.iterdir())} skills")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
