import json
from typing import Any, Dict, List, Optional


def get_variant_annotations(variant_id: str, json_path: str) -> Dict[str, Optional[str]]:
    """
    Load JSON from `json_path`, search all variants for `variant_id`,
    and return its key annotation scores/labels.

    Returns
    -------
    dict
        {
            "ClinVar":        <str | None>,
            "PHACTboost":     <str | None>,
            "AlphaMissense":  <str | None>,
            "GnomAD":         <str | None>,
        }

    Raises
    ------
    ValueError
        If the JSON structure is unexpected or the variant isn’t found.
    """
    # 1. Load the JSON file
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # 2. Gather all variants into one list
    variants: List[Dict[str, Any]] = []

    def _collect(rec: Dict[str, Any]) -> None:
        for key in ("variants", "Variants"):
            lst = rec.get(key)
            if isinstance(lst, list):
                variants.extend(lst)

    if isinstance(data, dict):
        _collect(data)
    elif isinstance(data, list):
        for record in data:
            if isinstance(record, dict):
                _collect(record)
    else:
        raise ValueError("JSON must be a dict or a list of dicts")

    # 3. Find the matching variant and extract annotations
    for var in variants:
        if var.get("VariantID") == variant_id:
            ann = var.get("annotations", {})
            # Safely dig through nested dicts
            clinvar = (
                ann.get("Variant databases", {}).get("ClinVar")
                if isinstance(ann.get("Variant databases"), dict)
                else None
            )
            path_pred = ann.get("Pathogenicity predictions", {})
            phact = path_pred.get("PHACTboost") if isinstance(path_pred, dict) else None
            alpha = path_pred.get("AlphaMissense") if isinstance(path_pred, dict) else None
            gnomad = (
                ann.get("Population allele frequency", {}).get(
                    "GnomAD joint allele frequency"
                )
                if isinstance(ann.get("Population allele frequency"), dict)
                else None
            )

            return {
                "ClinVar": clinvar,
                "PHACTboost": phact,
                "AlphaMissense": alpha,
                "GnomAD": gnomad,
            }

    # 4. Not found
    raise ValueError(f"VariantID '{variant_id}' not found")
result = get_variant_annotations("X-156001632-T-A", "mgs.usergenes_Wes_3687.json")
print(result)
# {'ClinVar': 'Benign', 'PHACTboost': 'None', 'AlphaMissense': 'N/A', 'GnomAD': '0.366891'}
