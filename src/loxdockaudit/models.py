"""Data models for LOXDockAudit."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ActiveSiteConfig:
    residues: list[dict]


@dataclass
class TargetResidue:
    chain: str
    resi: int
    resn: str
    atom: str = "NZ"


@dataclass
class DomainRange:
    start: int
    end: int
    name: str = ""


@dataclass
class ConstructConfig:
    construct_id: str
    construct_type: str
    models_dir: str
    receptor_chains: list[str]
    ligand_chains: list[str]
    active_site: ActiveSiteConfig
    target_residues: list[TargetResidue]
    lox_domain: DomainRange
    linker: DomainRange | None = None
    cbd_domain: DomainRange | None = None
    baseline_construct: str | None = None
    productive_distance_threshold: float = 8.0
    contact_distance_threshold: float = 4.0
    best_productive_rank_max: int = 3
    is_inactive: bool = False
    control_type: str | None = None


@dataclass
class PoseMetrics:
    rank: int
    model_path: str
    active_site_to_target_distance: float
    productive: bool
    lox_contacts: int
    cbd_contacts: int
    closest_active_site_resi: int | None = None
    warning: str | None = None


@dataclass
class ConstructSummary:
    construct_id: str
    productive_count: int
    total_poses: int
    best_distance: float
    best_productive_rank: int | None
    lox_contact_frequency: float
    cbd_contact_frequency: float
    cbd_coupling: float
    strict_pass: bool | None


@dataclass
class ControlEntry:
    construct_id: str
    control_type: str
    config_path: str


@dataclass
class StrictCriteria:
    pass_count_greater_than_baseline: bool = True
    best_distance_lte_baseline: bool = True
    best_productive_rank_max: int = 3
    cbd_must_contribute: bool = False
    must_beat_scrambled: bool = True
    must_beat_polyK: bool = True


@dataclass
class ScreenConfig:
    screen_id: str
    candidate_construct_id: str
    constructs: list[ControlEntry]
    strict_criteria: StrictCriteria


@dataclass
class ControlComparisonResult:
    construct_id: str
    control_type: str
    candidate_productive_count: int
    control_productive_count: int
    candidate_best_distance: float
    control_best_distance: float
    candidate_beats_control_count: bool
    candidate_beats_control_distance: bool
    conclusion: str


@dataclass
class ScreenResult:
    screen_id: str
    candidate_summary: ConstructSummary
    baseline_summary: ConstructSummary | None
    control_comparisons: list[ControlComparisonResult]
    criteria_results: dict[str, bool]
    overall_pass: bool
    decision_text: str
