# long-note-tex mode

Use this mode for durable technical documentation or a teaching-oriented scientific guide. Its signature is progressive explanation: establish the physical or conceptual picture before introducing inference, algorithms, and implementation details.

## Document architecture

Use only the sections that serve the subject, typically in this order:

1. Title, abstract, and table of contents.
2. Problem statement and desired map.
3. Physical or conceptual picture, with a labeled figure when useful.
4. Governing equations and term-by-term interpretation.
5. Known inputs, unknown quantities, and observation/data model.
6. Deterministic calibration or solution workflow before specialized probability or optimization language.
7. Progressive mathematical concepts and derivations.
8. Numerical implementation, algorithm, or pseudocode.
9. Identifiability, uncertainty, limitations, and validation.
10. Recommended implementation sequence, glossary, and optional appendices.

## Visual system

- `article`, 11 pt, approximately 0.9 inch margins.
- Blue primary accent plus restrained orange, green, and purple semantic accents.
- `microtype`, modest paragraph spacing, and no paragraph indentation.
- Use `\fcolorbox` callouts for one central idea or warning; do not turn every paragraph into a box.
- Use TikZ for explanatory process diagrams and physical schematics when prose alone is insufficient.
- Use `listings` for compact pseudocode with line numbers, wrapping, quiet background, and captions.
- Use `booktabs` tables and a glossary when the intended audience may not know the terminology.

## Explanatory sequence

- Begin with what the reader is trying to predict or decide.
- Give the physical meaning of every state and coefficient before the full equation.
- Read coupled equations term by term after displaying them.
- Explain deterministic calibration before prior/likelihood/posterior language.
- Separate coefficient space, state space, observation space, and numerical grid.
- Explain what uncertainty bands contain and what they omit.
- State what synthetic examples prove and what they do not prove.

## Equation quality

Define the forward problem completely: domain, state, parameters, forcing, initial conditions, boundary conditions, observation operator, and output. Use dimensional reference scales for positive transforms and dimensionless latent variables. Avoid bare symbolic objectives whose residual, norm, covariance, or regularization is undefined.

## Verification

- Compile twice and inspect the table of contents, figures, captions, algorithms, references, and appendices.
- Require zero fatal errors, undefined references, and overfull boxes.
- Check that floats appear near the introducing text and that callout boxes do not split awkwardly.
- For a very long note, confirm that repeated definitions were consolidated rather than copied across sections.

