#!/usr/bin/env python3
"""Build protected local HRS analysis datasets from an authorized RAND source.

This script is supplied without HRS data. It never modifies the source ZIP.
Run it only in a local non-AI environment using an authorized data product.
The derived outputs must remain protected and are excluded from this repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


MEMBER = "randhrs1992_2022v1.dta"
SOURCE_SHA256 = "354769dbf0975e5879971046dede3e05e21be4bce572541cc42cd9958ea425cc"
ITEM_STEMS = ("batha", "dressa", "eata", "beda", "walkra")
CONDITIONS = {
    "hypertension": ("hibp", "hibpe"),
    "diabetes": ("diab", "diabe"),
    "cancer": ("cancr", "cancre"),
    "lung_disease": ("lung", "lunge"),
    "heart_disease": ("heart", "hearte"),
    "stroke": ("strok", "stroke"),
    "psychiatric_problem": ("psych", "psyche"),
    "arthritis": ("arthr", "arthre"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def binary_or_nan(series: pd.Series) -> pd.Series:
    values = numeric(series)
    return values.where(values.isin([0, 1]))


def normalized_hhidpn(series: pd.Series) -> pd.Series:
    """Normalize HHIDPN only long enough to apply a local official data alert."""
    text = series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)
    return text.str.zfill(9)


def load_official_duplicate_drop_ids(path: Path) -> set[str]:
    """Load a local, authorized record-correction list without publishing identifiers.

    Record identifiers are deliberately absent from this code-only repository.
    An eligible user must create this one-identifier-per-line file locally from
    the applicable official HRS data alert and keep it outside version control.
    """
    values = {
        line.strip().zfill(9)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    if not values:
        raise RuntimeError("official duplicate-ID file is empty")
    if any(not re.fullmatch(r"\d{9}", value) for value in values):
        raise RuntimeError("official duplicate-ID file must contain one 9-digit identifier per line")
    return values


def cutoff_state(ever: pd.Series, flag: pd.Series) -> pd.Series:
    ever_num = binary_or_nan(ever)
    flag_num = numeric(flag)
    result = ever_num.copy()
    result.loc[flag_num.eq(6)] = 1.0
    result.loc[flag_num.isin([4, 5])] = np.nan
    return result


def build_pair(
    source: pd.DataFrame,
    baseline: int,
    followup: int,
    source_path: Path,
    official_duplicate_drop_ids: set[str],
) -> tuple[pd.DataFrame, dict]:
    out = pd.DataFrame({"source_row": np.arange(len(source), dtype=np.int64)})
    out["official_duplicate_drop"] = normalized_hhidpn(source["hhidpn"]).isin(
        official_duplicate_drop_ids
    ).astype(int)

    def take(field: str) -> pd.Series:
        return numeric(source[field])

    out["base_iwstat"] = take(f"r{baseline}iwstat")
    out["follow_iwstat"] = take(f"r{followup}iwstat")
    out["age"] = take(f"r{baseline}agey_b")
    out["sex_code"] = take("ragender")
    out["hispanic_code"] = take("rahispan")
    out["race_code"] = take("raracem")
    out["education_years"] = take("raedyrs")
    out["marital_code"] = take(f"r{baseline}mstat")
    out["base_proxy"] = take(f"r{baseline}proxy")
    out["follow_proxy"] = take(f"r{followup}proxy")
    out["follow_nursing_home"] = take(f"r{followup}nhmliv")
    out["bmi"] = take(f"r{baseline}bmi")
    out["base_weight"] = take(f"r{baseline}wtcrnh")
    out["raestrat"] = take("raestrat")
    out["raehsamp"] = take("raehsamp")
    mode_field = f"r{followup}cogmode"
    out["follow_cogmode"] = take(mode_field) if mode_field in source else np.nan

    baseline_items = []
    follow_items = []
    for stem in ITEM_STEMS:
        b = binary_or_nan(source[f"r{baseline}{stem}"])
        f = binary_or_nan(source[f"r{followup}{stem}"])
        out[f"base_{stem}"] = b
        out[f"follow_{stem}"] = f
        baseline_items.append(b)
        follow_items.append(f)
    baseline_matrix = pd.concat(baseline_items, axis=1)
    follow_matrix = pd.concat(follow_items, axis=1)
    baseline_complete = baseline_matrix.notna().all(axis=1)
    follow_complete = follow_matrix.notna().all(axis=1)
    out["eligible_base"] = (
        out["age"].ge(65)
        & baseline_complete
        & baseline_matrix.eq(0).all(axis=1)
        & out["base_weight"].gt(0)
        & out["base_iwstat"].eq(1)
        & out["official_duplicate_drop"].eq(0)
    ).astype(int)
    out["full_weighted_frame"] = (
        out["base_iwstat"].eq(1)
        & out["base_weight"].gt(0)
        & out["raestrat"].notna()
        & out["raehsamp"].notna()
        & out["official_duplicate_drop"].eq(0)
    ).astype(int)

    any_positive = follow_matrix.eq(1).any(axis=1)
    known_negative = follow_complete & follow_matrix.eq(0).all(axis=1)
    y_func = pd.Series(np.nan, index=source.index, dtype=float)
    y_func.loc[any_positive] = 1.0
    y_func.loc[known_negative] = 0.0
    out["y_function"] = y_func
    out["y_function_complete_items"] = np.where(
        follow_complete, follow_matrix.eq(1).any(axis=1).astype(float), np.nan
    )
    out["follow_alive"] = out["follow_iwstat"].isin([1, 4]).astype(int)
    out["follow_death"] = out["follow_iwstat"].eq(5).astype(int)

    final_cols = []
    cutoff_cols = []
    ambiguous_cols = []
    other_unknown_cols = []
    condition_qc = {}
    for condition, (stem, ever_stem) in CONDITIONS.items():
        ever = binary_or_nan(source[f"r{baseline}{ever_stem}"])
        flag = take(f"r{baseline}{stem}f")
        raw = take(f"r{baseline}{stem}")
        question = take(f"r{baseline}{stem}q")
        cutoff = cutoff_state(ever, flag)
        ambiguous = flag.isin([4, 5])
        other_unknown = cutoff.isna() & ~ambiguous
        out[f"final_{condition}"] = ever
        out[f"cutoff_{condition}"] = cutoff
        out[f"flag_{condition}"] = flag
        out[f"raw_{condition}"] = raw
        out[f"question_{condition}"] = question
        final_cols.append(f"final_{condition}")
        cutoff_cols.append(f"cutoff_{condition}")
        ambiguous_cols.append(ambiguous)
        other_unknown_cols.append(other_unknown)
        condition_qc[condition] = {
            "flag4": int(flag.eq(4).sum()),
            "flag5": int(flag.eq(5).sum()),
            "flag6": int(flag.eq(6).sum()),
        }

    final_matrix = out[final_cols]
    cutoff_matrix = out[cutoff_cols]
    ambiguous_matrix = pd.concat(ambiguous_cols, axis=1)
    other_unknown_matrix = pd.concat(other_unknown_cols, axis=1)
    out["final_exposure_known"] = final_matrix.notna().all(axis=1).astype(int)
    out["cutoff_exposure_known"] = cutoff_matrix.notna().all(axis=1).astype(int)
    out["exposure_identified"] = (
        final_matrix.notna().all(axis=1) & cutoff_matrix.notna().all(axis=1)
    ).astype(int)
    out["ambiguous_condition_count"] = ambiguous_matrix.sum(axis=1).astype(int)
    out["other_unknown_condition_count"] = other_unknown_matrix.sum(axis=1).astype(int)
    out["extreme_scenario_eligible"] = (
        final_matrix.notna().all(axis=1) & ~other_unknown_matrix.any(axis=1)
    ).astype(int)
    out["final_count"] = final_matrix.sum(axis=1, min_count=len(CONDITIONS))
    out["cutoff_count"] = cutoff_matrix.sum(axis=1, min_count=len(CONDITIONS))
    cutoff_zero = cutoff_matrix.fillna(0)
    cutoff_one = cutoff_matrix.fillna(1)
    out["cutoff_count_all0"] = cutoff_zero.sum(axis=1)
    out["cutoff_count_all1"] = cutoff_one.sum(axis=1)

    out["primary_survivor_target"] = (
        out["eligible_base"].eq(1)
        & out["exposure_identified"].eq(1)
        & out["follow_alive"].eq(1)
    ).astype(int)
    out["primary_outcome_observed"] = (
        out["primary_survivor_target"].eq(1) & out["y_function"].notna()
    ).astype(int)
    out["adverse_state"] = np.where(
        out["follow_death"].eq(1),
        1.0,
        np.where(out["follow_alive"].eq(1) & out["y_function"].notna(), out["y_function"], np.nan),
    )
    out["adverse_target"] = (
        out["eligible_base"].eq(1)
        & out["exposure_identified"].eq(1)
        & out["follow_iwstat"].isin([1, 4, 5])
    ).astype(int)

    eligible = out["eligible_base"].eq(1)
    identified = eligible & out["exposure_identified"].eq(1)
    survivors = identified & out["follow_alive"].eq(1)
    observed = survivors & out["y_function"].notna()
    changed = identified & out["final_count"].ne(out["cutoff_count"])
    qc = {
        "source": str(source_path),
        "source_sha256_expected": SOURCE_SHA256,
        "baseline_wave": baseline,
        "followup_wave": followup,
        "source_rows": int(len(out)),
        "official_duplicate_alert_rows_excluded": int(out["official_duplicate_drop"].sum()),
        "full_weighted_frame": int(out["full_weighted_frame"].sum()),
        "eligible_base": int(eligible.sum()),
        "eligible_exposure_identified": int(identified.sum()),
        "eligible_with_any_f4_or_f5": int((eligible & out["ambiguous_condition_count"].gt(0)).sum()),
        "eligible_other_exposure_unknown": int((eligible & out["other_unknown_condition_count"].gt(0)).sum()),
        "identified_followup_survivors": int(survivors.sum()),
        "identified_survivors_outcome_observed": int(observed.sum()),
        "identified_survivors_outcome_positive": int((observed & out["y_function"].eq(1)).sum()),
        "identified_survivors_outcome_negative": int((observed & out["y_function"].eq(0)).sum()),
        "identified_policy_count_changed_baseline": int(changed.sum()),
        "identified_policy_count_change_distribution": {
            str(int(key)): int(value)
            for key, value in (
                out.loc[identified, "cutoff_count"] - out.loc[identified, "final_count"]
            ).value_counts().sort_index().items()
        },
        "deaths_in_identified_baseline": int((identified & out["follow_death"].eq(1)).sum()),
        "alive_unknown_in_primary_target": int((survivors & out["y_function"].isna()).sum()),
        "condition_flag_counts_all_rows": condition_qc,
    }
    return out, qc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("derived"))
    parser.add_argument(
        "--source",
        type=Path,
        required=True,
        help="Local path to the authorized RAND HRS 1992–2022 v1 STATA ZIP archive",
    )
    parser.add_argument(
        "--official-duplicate-id-file",
        type=Path,
        required=True,
        help=(
            "Local, untracked one-identifier-per-line correction list created from "
            "the applicable official HRS data alert"
        ),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_path = args.source.expanduser().resolve()
    official_duplicate_drop_ids = load_official_duplicate_drop_ids(
        args.official_duplicate_id_file.expanduser().resolve()
    )

    fields = {
        "hhidpn", "ragender", "rahispan", "raracem", "raedyrs", "raestrat", "raehsamp"
    }
    for baseline, followup in ((14, 15), (15, 16)):
        fields.update(
            {
                f"r{baseline}iwstat", f"r{followup}iwstat",
                f"r{baseline}agey_b", f"r{baseline}proxy", f"r{followup}proxy",
                f"r{followup}nhmliv", f"r{baseline}mstat", f"r{baseline}bmi",
                f"r{baseline}wtcrnh",
            }
        )
        for stem in ITEM_STEMS:
            fields.add(f"r{baseline}{stem}")
            fields.add(f"r{followup}{stem}")
        for stem, ever_stem in CONDITIONS.values():
            fields.update(
                {
                    f"r{baseline}{stem}", f"r{baseline}{stem}q",
                    f"r{baseline}{stem}f", f"r{baseline}{ever_stem}",
                }
            )
    fields.add("r15cogmode")

    chunks = []
    with zipfile.ZipFile(source_path) as archive:
        with archive.open(MEMBER) as stream:
            for chunk in pd.read_stata(
                stream,
                iterator=True,
                columns=sorted(fields),
                chunksize=5000,
                convert_categoricals=False,
            ):
                chunks.append(chunk)
    source = pd.concat(chunks, ignore_index=True)

    manifest = {
        "source": str(source_path),
        "source_sha256_expected": SOURCE_SHA256,
        "source_sha256_observed": sha256(source_path),
        "member": MEMBER,
        "source_rows": int(len(source)),
        "raw_identifier_retained": False,
        "official_duplicate_drop_ids_count": len(official_duplicate_drop_ids),
        "outputs": [],
    }
    if manifest["source_sha256_observed"] != SOURCE_SHA256:
        raise RuntimeError("source ZIP hash differs from the Q1 lock")

    for baseline, followup in ((14, 15), (15, 16)):
        data, qc = build_pair(
            source,
            baseline,
            followup,
            source_path,
            official_duplicate_drop_ids,
        )
        data_path = args.output_dir / f"q2_hrs_w{baseline}_w{followup}.csv.gz"
        qc_path = args.output_dir / f"q2_hrs_w{baseline}_w{followup}_qc.json"
        data.to_csv(data_path, index=False, compression="gzip")
        qc_path.write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        manifest["outputs"].append(
            {
                "path": str(data_path),
                "sha256": sha256(data_path),
                "qc_path": str(qc_path),
                "qc_sha256": sha256(qc_path),
            }
        )

    manifest_path = args.output_dir / "q2_derived_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
