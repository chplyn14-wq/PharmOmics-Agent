"""Tests for pharmomics.cli.analyze — build_experiment_design."""

from __future__ import annotations

import pytest

from pharmomics.cli.analyze import build_experiment_design
from pharmomics.experiment.enums import GroupRole


class TestBuildExperimentDesignBasic:
    """Verify metadata → ExperimentDesign construction."""

    def _sample_ids(self) -> list[str]:
        return ["ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"]

    def _conditions(self) -> dict[str, str]:
        return {
            "ctrl_1": "DMSO",
            "ctrl_2": "DMSO",
            "ctrl_3": "DMSO",
            "trt_1": "osi_DTP",
            "trt_2": "osi_DTP",
            "trt_3": "osi_DTP",
        }

    def test_two_groups_created(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        assert len(design.groups) == 2

    def test_group_ids_are_slugified(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        group_ids = {g.group_id for g in design.groups}
        assert "dmso" in group_ids
        assert "osi_dtp" in group_ids

    def test_group_labels_preserve_condition(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        labels = {g.label for g in design.groups}
        assert "DMSO" in labels
        assert "osi_DTP" in labels

    def test_samples_assigned_to_correct_groups(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        for sample in design.samples:
            if sample.sample_id.startswith("ctrl"):
                assert sample.group_id == "dmso"
            else:
                assert sample.group_id == "osi_dtp"

    def test_control_group_has_control_role(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        ctrl_group = [g for g in design.groups if g.group_id == "dmso"][0]
        assert ctrl_group.role == GroupRole.CONTROL

    def test_treatment_group_has_treatment_role(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        trt_group = [g for g in design.groups if g.group_id == "osi_dtp"][0]
        assert trt_group.role == GroupRole.TREATMENT

    def test_single_contrast_created(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        assert len(design.contrasts) == 1
        c = design.contrasts[0]
        assert c.comparison_group_id == "osi_dtp"
        assert c.reference_group_id == "dmso"

    def test_contrast_id_is_slugified(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        assert design.contrasts[0].contrast_id == "osi_dtp_vs_dmso"

    def test_factor_condition_created(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        assert len(design.factors) == 1
        assert design.factors[0].factor_id == "condition"
        assert set(design.factors[0].levels or []) == {"DMSO", "osi_DTP"}

    def test_sample_order_preserved(self) -> None:
        sids = self._sample_ids()
        design = build_experiment_design(
            sample_ids=sids,
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        assert [s.sample_id for s in design.samples] == sids

    def test_factor_values_set(self) -> None:
        design = build_experiment_design(
            sample_ids=self._sample_ids(),
            conditions=self._conditions(),
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        for sample in design.samples:
            assert (
                sample.factor_values["condition"]
                == self._conditions()[sample.sample_id]
            )

    def test_observational_group_role(self) -> None:
        conditions = {
            "ctrl_1": "DMSO",
            "trt_1": "osi_DTP",
            "other_1": "other_cond",
        }
        design = build_experiment_design(
            sample_ids=["ctrl_1", "trt_1", "other_1"],
            conditions=conditions,
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        other_group = [g for g in design.groups if g.group_id == "other_cond"][0]
        assert other_group.role == GroupRole.OBSERVATIONAL


class TestBuildExperimentDesignErrors:
    """Verify error handling for invalid input."""

    def _sample_ids(self) -> list[str]:
        return ["ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"]

    def _conditions(self) -> dict[str, str]:
        return {
            "ctrl_1": "DMSO",
            "ctrl_2": "DMSO",
            "ctrl_3": "DMSO",
            "trt_1": "osi_DTP",
            "trt_2": "osi_DTP",
            "trt_3": "osi_DTP",
        }

    def test_missing_condition_raises(self) -> None:
        with pytest.raises(ValueError, match="no condition"):
            build_experiment_design(
                sample_ids=["S1", "S2"],
                conditions={"S1": "ctrl"},
                contrast_control="ctrl",
                contrast_treatment="trt",
            )

    def test_control_has_no_samples_raises(self) -> None:
        conditions = {"S1": "trt", "S2": "trt"}
        with pytest.raises(ValueError, match="Control condition"):
            build_experiment_design(
                sample_ids=["S1", "S2"],
                conditions=conditions,
                contrast_control="ctrl",
                contrast_treatment="trt",
            )

    def test_treatment_has_no_samples_raises(self) -> None:
        conditions = {"S1": "ctrl", "S2": "ctrl"}
        with pytest.raises(ValueError, match="Treatment condition"):
            build_experiment_design(
                sample_ids=["S1", "S2"],
                conditions=conditions,
                contrast_control="ctrl",
                contrast_treatment="trt",
            )

    def test_empty_sample_ids_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            build_experiment_design(
                sample_ids=[],
                conditions={},
                contrast_control="ctrl",
                contrast_treatment="trt",
            )

    def test_extra_condition_in_dict_ignored(self) -> None:
        # Extra entries in conditions dict that aren't in sample_ids
        # don't affect the design (cross-validation is done upstream)
        design = build_experiment_design(
            sample_ids=["S1", "S2"],
            conditions={"S1": "ctrl", "S2": "trt", "S3": "other"},
            contrast_control="ctrl",
            contrast_treatment="trt",
        )
        assert len(design.samples) == 2
        group_ids = {g.group_id for g in design.groups}
        # "other" condition still creates a group since it's in conditions.values()
        assert "other" in group_ids
