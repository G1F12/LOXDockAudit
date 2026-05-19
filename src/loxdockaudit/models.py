"""Data models for LOXDockAudit."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


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
class QCConfig:
    his_resi: list[int]
    lys_resi: int
    tyr_resi: int
    disulfide_pairs: list[tuple[int, int]] = field(default_factory=list)
    alphafold_pdb_path: str | None = None
    pae_json_path: str | None = None
    domain_a_resi: list[int] | None = None
    domain_b_resi: list[int] | None = None
    plddt_threshold: float = 70.0
    his_max_ca_distance: float = 10.0
    lys_tyr_max_cb_distance: float = 12.0
    disulfide_max_sg_distance: float = 2.5


# v0.4: construct-level orientation scoring configuration.
@dataclass
class OrientationConfig:
    enabled: bool = False
    threshold_deg: float = 90.0
    sidechain_atoms: list[str] = field(default_factory=lambda: ["CB", "NZ"])
    activesite_centroid_atoms: list[str] = field(default_factory=lambda: ["CA"])


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
    structural_qc: QCConfig | None = None
    orientation: OrientationConfig = field(default_factory=OrientationConfig)  # v0.4


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
    orientation_angle_deg: float | None = None  # v0.4
    orientation_productive: bool | None = None  # v0.4
    fully_productive: bool | None = None  # v0.4
    orientation_enabled: bool = False  # v0.4


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
    orientation_enabled: bool = False  # v0.4
    fully_productive_count: int = 0  # v0.4
    best_orientation_angle_deg: float | None = None  # v0.4
    best_fully_productive_rank: int | None = None  # v0.4


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
    # v0.4: None means the optional orientation criterion was absent from YAML.
    fully_productive_count_greater_than_baseline: bool | None = None
    best_orientation_angle_max_deg: float | None = None
    best_fully_productive_rank_max: int | None = None


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
    criteria_results: dict[str, Any]
    overall_pass: bool
    decision_text: str


@dataclass
class StructuralQCResult:
    construct_id: str
    his_triad_pass: bool
    lys_tyr_pass: bool
    disulfide_pass: bool
    active_site_accessible: bool
    plddt_pass: bool | None
    pae_pass: bool | None
    fold_qc_pass: bool
    fold_corrupted: bool
    warnings: list[str]
    details: dict
