"""PharmOmics analysis specification domain.

AnalysisSpecification represents analytical intent — what question to ask
of the data — separate from how the experiment was structured or what was
measured.

See ADR 0007 for the rationale: contrasts belong to AnalysisSpecification,
not ExperimentDesign.
"""

from pharmomics.analysis.contrast import ResolvedContrast, resolve_contrast
from pharmomics.analysis.matrix_slice import (
    DifferentialInput,
    MatrixSlice,
    MatrixSliceError,
    prepare_differential_inputs,
)
from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.analysis.validation import validate_analysis_specification

__all__ = [
    "AnalysisSpecification",
    "DifferentialInput",
    "MatrixSlice",
    "MatrixSliceError",
    "ResolvedContrast",
    "prepare_differential_inputs",
    "resolve_contrast",
    "validate_analysis_specification",
]
