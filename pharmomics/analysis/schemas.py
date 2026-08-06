from dataclasses import dataclass, field


@dataclass(frozen=True)
class AnalysisSpecification:
    analysis_type: str
    factor_references: list[str] = field(default_factory=list)
    contrast_references: list[str] = field(default_factory=list)
    parameters: dict[str, object] = field(default_factory=dict)
