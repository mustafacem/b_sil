import json
from pathlib import Path
from typing import Dict, List, Any


def _flatten_inheritance(val: Any) -> List[str]:
    """Helper – turn whatever is stored under an “Inheritance” field into a flat list of strings."""
    if val is None:
        return []
    out: List[str] = []
    if isinstance(val, str):
        out.append(val)
    elif isinstance(val, list):
        out.extend(x for x in val if x)
    elif isinstance(val, dict):
        for v in val.values():
            out.extend(_flatten_inheritance(v))
    return out


def get_inheritance_for_variant(
    json_path: str | Path,
    variant_id: str,
) -> Dict[str, List[str]]:
    """
    Return `{phenotype_name_en: [inheritance, …], …}` for *variant_id*.

    ‣ Name comes from `phenotype["Name"]["en"]`  
      – falls back to the Turkish name, then to the raw phenotypeID.  
    ‣ If the variant does not list phenotype objects, the function tries
      the lightweight `Phenotype` + `PhenotypeName` arrays.

    Raises
    ------
    FileNotFoundError | json.JSONDecodeError | ValueError
    """
    json_path = Path(json_path)

    with json_path.open("r", encoding="utf-8") as fh:
        data: Any = json.load(fh)
    if isinstance(data, dict):
        data = [data]

    for gene_block in data:
        for variant in gene_block.get("Variants", []):
            if variant.get("VariantID") != variant_id:
                continue

            result: Dict[str, List[str]] = {}

            # ── Preferred: rich phenotype objects ───────────────────────────
            for ph in variant.get("phenotypes", []):
                pid = ph.get("phenotypeID") or ph.get("_id", {}).get("$oid") or ""
                name_en = (ph.get("Name", {}) or {}).get("en", "") or \
                          (ph.get("Name", {}) or {}).get("tr", "") or pid
                result[name_en] = _flatten_inheritance(ph.get("Inheritance"))

            # ── Fallback: plain arrays (Phenotype / PhenotypeName) ──────────
            if not result and "Phenotype" in variant:
                ids: List[str] = variant["Phenotype"]
                names: List[str] = variant.get("PhenotypeName", [])
                inh:  List[Any]  = variant.get("Inheritance", [])

                for idx, pid in enumerate(ids):
                    raw_name = names[idx] if idx < len(names) else ""
                    display_name = raw_name if raw_name and raw_name.lower() != "none" else pid
                    inherit_list = _flatten_inheritance(inh[idx] if idx < len(inh) else None)
                    result[display_name] = inherit_list

            return result

    raise ValueError(f"VariantID '{variant_id}' not found in {json_path}")


# ─────────────────────────── Example ────────────────────────────
if __name__ == "__main__":
    path = "mgs.usergenes_Wes_3687.json"
    vid  = "X-154905041-CA-C"

    pheno_inheritance = get_inheritance_for_variant(path, vid)

    # If you only need the first inheritance term for each phenotype:
    simple = {name: vals[0] if vals else "" for name, vals in pheno_inheritance.items()}
    print(simple)
    items_str = [f"{value}: {key}" for key, value in simple.items()]
    output_str = ", ".join(items_str)
    print(output_str.replace(":", " for "))