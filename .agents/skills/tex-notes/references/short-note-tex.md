# short-note-tex mode

Use this mode for a dense but readable decision brief. Its signature is one overview page followed by self-contained method pages, with no method spilling onto a second page.

## Page architecture

1. **Overview page:** title, date/status, objective, common reference model, notation, a compact route/method table, and the decision order.
2. **One page per method:** method title and role, purpose, forward equations, inputs/outputs, implementation requirements, and merit/risk. Begin every method with `\clearpage` and `\MethodTitle{...}{...}`.
3. **Optional final page:** comparison table, shared evidence gates, recommended sequence, and selected references.

If a method page overflows, shorten prose, move shared definitions to the overview, combine repetitive bullets, and simplify displays. Do not solve overflow by making body text uncomfortably small.

## Visual system

- `article`, 10 pt, letter paper, approximately 0.62 inch margins.
- Dark blue method titles, quiet gray role labels, plain page numbering.
- Compact paragraph and display spacing; no decorative cover page.
- `booktabs` and ragged-right `tabularx` columns for readable narrow tables.
- Box only the single equation or conclusion that deserves primary emphasis.
- Keep equation labels semantic and cross-references resolved.

## Scientific writing pattern

For each method, answer in this order:

1. What problem does it solve?
2. What is its state and forward equation?
3. Which quantities are known inputs, learned coefficients, forcing, initial data, and outputs?
4. What identities or invariants should the formulation preserve?
5. What must be implemented or measured?
6. What is the realistic merit, failure mode, and evidence threshold?

Use one notation across every page. When a symbol could conflict across methods, rename it once in the overview rather than redefining it locally.

## Verification

- Compile twice in a temporary directory.
- Confirm every method begins on a new page and occupies no more than one compiled page.
- Require zero fatal errors, undefined references, and overfull boxes. Underfull warnings in narrow comparison tables should be removed with ragged-right column types rather than ignored when they visibly harm typography.

