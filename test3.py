import json
from typing import List, Union, Any, Dict, Set


# ───────────────────────────── helpers ────────────────────────────────────
def _extract_records(obj: Any) -> List[dict]:
    if isinstance(obj, list) and all(isinstance(x, dict) for x in obj):
        return obj
    if isinstance(obj, dict):
        for v in obj.values():
            if isinstance(v, list) and v and isinstance(v[0], dict):
                return v
    return [obj] if isinstance(obj, dict) else []


def _find_descriptions(obj: Any):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k.lower() == "description":
                yield v
            else:
                yield from _find_descriptions(v)
    elif isinstance(obj, list):
        for item in obj:
            yield from _find_descriptions(item)






def _has_pathogenic_label(variant: Dict[str, Any]) -> bool:
    top = variant.get("Final_label", [])
    if isinstance(top, str):
        top = [top]
    if any("pathogenic" in str(lbl).lower() for lbl in top):
        return True

    for item in variant.get("pheno_inh", []):
        lbl = item.get("Final_label")
        if isinstance(lbl, str) and "pathogenic" in lbl.lower():
            return True
        if isinstance(lbl, list) and any("pathogenic" in str(x).lower() for x in lbl):
            return True
    return False





def extract_en_pathogenic_names(rec: Dict[str, Any]) -> List[str]:
    """
    Only return Name['en'] for phenotypes that belong to variants
    whose Final_label mentions 'pathogenic'.
    """
    seen: Set[str] = set()
    for var in rec.get("Variants", []):
        if not _has_pathogenic_label(var):           # ← NEW guard
            continue
        for pheno in var.get("phenotypes", []):
            name_en = pheno.get("Name", {}).get("en")
            if name_en:
                seen.add(name_en.strip())
    return sorted(seen)


def filter_gene_names_by_description(
    data: Union[str, List[dict]],
    significance_threshold: float = 0.8,
    keywords: List[str] = ["neurological", "neuron", "brain"]
) -> List[str]:
    """
    A record is kept only if:
      • record["significance"] ≥ threshold
      • some Description text contains a keyword
      • it contains at least one *pathogenic* variant (checked implicitly
        because extract_en_pathogenic_names() will return [] otherwise)
    Returned list consists solely of phenotype Name['en'] strings coming
    from pathogenic variants.
    """
    # load JSON if a filename was supplied
    if isinstance(data, str):
        with open(data, "r", encoding="utf-8") as f:
            raw = json.load(f)
    else:
        raw = data

    records = _extract_records(raw)
    kws = {kw.lower() for kw in keywords}
    result: List[str] = []

    for rec in records:
        if float(rec.get("significance", 0)) < significance_threshold:
            continue

        if not any(
            any(kw in ((" ".join(str(v) for v in desc.values())
                        if isinstance(desc, dict) else str(desc)).lower())
                for kw in kws)
            for desc in _find_descriptions(rec)
        ):
            continue

        names = extract_en_pathogenic_names(rec)     # ← uses new function
        if names:
            print("-----------------------------------------------------------------")
            print(names)
            result.extend(names)

    return result


# quick test
if __name__ == "__main__":
    genes = filter_gene_names_by_description("mgs.usergenes_Wes_3687.json")
    print(", ".join(genes))