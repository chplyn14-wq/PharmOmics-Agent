\# PharmOmics Agent



\## Project goal



Build a human-in-the-loop research assistant for drug-response

transcriptomics. The system performs deterministic QC and statistical

analysis, retrieves PubMed evidence, generates evidence-grounded

hypotheses, and proposes high-level validation experiments.



\## Scientific rules



\- Never fabricate genes, statistics, citations, PMIDs or experimental results.

\- Numerical results must come from deterministic analysis code.

\- LLMs may interpret results but must never modify calculated values.

\- Every literature-supported claim must include a valid PMID.

\- Distinguish supporting, contradictory and indirect evidence.

\- Every hypothesis must include limitations and falsification criteria.

\- Do not describe outputs as clinically validated.

\- Require human approval before producing experiment recommendations.



\## Engineering rules



\- Use typed Python and Pydantic models.

\- Write tests before or alongside each feature.

\- Keep analysis tools separate from LLM agents.

\- Store every run in a versioned run directory.

\- Do not commit API keys, large raw data or generated run files.

\- Every public function must have a docstring.

\- Do not silently catch errors.

\- Prefer small, reviewable commits.



\## Commands



\- Run tests: pytest

\- Run CLI: uv run pharmomics --help

\- Run formatting: ruff format .

\- Run linting: ruff check .

