---
name: paper-style-rewrite
description: Rewrite, diagnose, and tighten scientific paper prose with strict sentence-level gates for concision, technical precision, claim alignment, and non-AI-polished professor-like voice. Use when Codex is asked to rewrite a paragraph, section, abstract, introduction, result, discussion, conclusion, figure caption, reviewer response, or paper draft text; when the user asks for journal-ready scientific prose; when the user mentions avoiding AI-generated style; or when Gomez Lab paper-writing rules, claim alignment, no-waste sentences, or technical tightness should govern the output.
---

# Paper Style Rewrite

Use this skill to rewrite scientific paper prose without making it sound generically polished or AI-generated. The goal is sober, exact, professor-like writing in which every sentence contributes a distinct technical idea.

## Workflow

1. Read project-local context when available: `AGENTS.md`, `memory.md`, `rulebook.md`, `harness.md`, and the active draft context named there.
2. If project-local files are absent, read `references/rewrite-gates.md`.
3. Identify the passage's technical job before rewriting.
4. Rewrite without adding unsupported scientific claims.
5. Apply the no-waste, technical tightness, and professor-voice gates before returning text.
6. Return the rewritten passage first. Do not show the audit unless requested.

## Required Behavior

- Preserve the user's scientific meaning unless conceptual revision is requested.
- Do not edit files unless the user explicitly asks.
- Avoid first-person paper prose such as "we next ask", "we show", or "we test" unless requested.
- Remove filler recaps, announcement sentences, repeated ideas, and decorative transitions.
- Replace vague qualifiers with quantities, named baselines, explicit comparators, or mechanisms.
- Keep LaTeX notation, citations, and figure references intact unless the user asks to change them.
- Avoid promotional or stock phrases such as "paves the way", "underscores", "offers valuable insights", and "represents a significant advancement".

## Output Modes

For a simple rewrite request, return only the rewritten text unless a short note is needed to flag ambiguity.

For a rewrite plus diagnosis request, return the rewritten text first, then a concise list of issues fixed.

For a framing or logic question, answer the conceptual issue directly and suggest replacement wording only where useful.
