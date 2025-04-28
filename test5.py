import json
from pathlib import Path
from typing import Dict, List, Any


def _flatten_inheritance(val: Any) -> List[str]:
    """
    Helper – turn whatever is stored under an “Inheritance”
    field into a flat list of non-empty strings.
    """
    if val is None:
        return []

    # The gene/phenotype blocks sometimes store
    #  ▸ {"en": ["Autosomal dominant"], "tr": ["Otozomal dominant"]}
    #  ▸ ["Autosomal dominant"]
    #  ▸ "Autosomal dominant"
    out: List[str] = []
    if isinstance(val, str):
        out.append(val)
    elif isinstance(val, list):
        out.extend(x for x in val if x)
    elif isinstance(val, dict):
        for v in val.values():
            out.extend(_flatten_inheritance(v))
    return out


def get_inheritance_for_variant(json_path: str | Path, variant_id: str) -> Dict[str, List[str]]:
    """
    Return the inheritance mode(s) of *each* phenotype linked to a given VariantID.

    Parameters
    ----------
    json_path : str | Path
        Path to the JSON file on disk.
    variant_id : str
        The `VariantID` you are looking for (e.g. ``"X-154905041-CA-C"``).

    Returns
    -------
    dict
        ``{phenotypeID: ["Inheritance 1", "Inheritance 2", …]}``

        ▸ If the variant is annotated with phenotypes but *none* of them
          specify inheritance, the dict will map each phenotypeID to
          ``[]``.  
        ▸ If the variant itself has no phenotype annotations, an empty
          dict is returned.

    Raises
    ------
    FileNotFoundError
        *json_path* does not exist.
    json.JSONDecodeError
        File is not valid JSON.
    ValueError
        *variant_id* was not found anywhere in the file.
    """
    json_path = Path(json_path)

    # 1. Load & normalise the top-level JSON
    with json_path.open("r", encoding="utf-8") as fh:
        data: Any = json.load(fh)

    # The file is a list of gene blocks (as in your example).  If a single
    # object was given instead, wrap it so the same code works.
    if isinstance(data, dict):
        data = [data]

    # 2. Walk every gene-block → every variant until we hit the
    #    VariantID we’re after
    for gene_block in data:
        for variant in gene_block.get("Variants", []):
            if variant.get("VariantID") != variant_id:
                continue

            # 3. Build the result for this variant
            result: Dict[str, List[str]] = {}

            # Preferred source: the rich `phenotypes` objects the pipeline
            # attaches to some variants
            for ph in variant.get("phenotypes", []):
                pid = ph.get("phenotypeID") or ph.get("_id", {}).get("$oid")
                result[pid] = _flatten_inheritance(ph.get("Inheritance"))

            # Fallback path – some variants list phenotype IDs only and
            # put per-phenotype inheritance inside `variant["Inheritance"]`
            # at the same index position.
            if not result and "Phenotype" in variant:
                ids: List[str] = variant["Phenotype"]
                inh: List[Any] = variant.get("Inheritance", [])
                for idx, pid in enumerate(ids):
                    result[pid] = _flatten_inheritance(inh[idx] if idx < len(inh) else None)

            return result  # we’re done – stop searching

    # 4. Variant not found anywhere
    raise ValueError(f"VariantID '{variant_id}' not found in {json_path}")


# ───────────────────────────── Example ─────────────────────────────
if __name__ == "__main__":
    path = "mgs.usergenes_Wes_3687.json"          # the file you showed
    vid  = "X-154905041-CA-C"       
    dict1 = get_inheritance_for_variant(path, vid)
    result_dict1 = {key: value[0] for key, value in dict1.items()}              # pick any VariantID present
    print(result_dict1)

