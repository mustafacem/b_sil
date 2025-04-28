import json
from typing import List, Union, Any


def _extract_records(obj: Any) -> List[dict]:
    # if it’s already a list of dicts, use it; otherwise look for a single list‐of‐dicts value
    if isinstance(obj, list) and all(isinstance(x, dict) for x in obj):
        return obj
    if isinstance(obj, dict):
        for v in obj.values():
            if isinstance(v, list) and v and isinstance(v[0], dict):
                return v
    return [obj] if isinstance(obj, dict) else []

def _find_descriptions(obj: Any):
    # yield every value under any “Description” key (case‐insensitive)
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k.lower() == "description":
                yield v
            else:
                yield from _find_descriptions(v)
    elif isinstance(obj, list):
        for item in obj:
            yield from _find_descriptions(item)

from typing import Dict, List, Set, Any

def extract_en_names(rec: Dict[str, Any]) -> List[str]:
    """
    Return a deduplicated list of phenotype names (`Name["en"]`)
    contained anywhere inside `rec["Variants"][*]["phenotypes"]`.
    """
    seen: Set[str] = set()

    # walk through every variant → every phenotype
    for var in rec.get("Variants", []):
        for pheno in var.get("phenotypes", []):
            en_name = pheno.get("Name", {}).get("en")
            if en_name:                       # skip empty or None
                seen.add(en_name.strip())

    return sorted(seen)                      # nice, predictable order



def filter_gene_names_by_description(
    data: Union[str, List[dict]],
    significance_threshold: float = 0.8,
    keywords: List[str] = ["neurological", "neuron", "brain"]
) -> List[str]:
    """
    Returns just the geneProperties.gene_name of records where:
      - record['significance'] >= significance_threshold
      - any Description text contains one of the keywords
    """
    # load JSON if a filename was passed
    if isinstance(data, str):
        with open(data, 'r', encoding='utf-8') as f:
            raw = json.load(f)
    else:
        raw = data

    records = _extract_records(raw)
    kws = {kw.lower() for kw in keywords}
    result = []

    for rec in records:
        if float(rec.get("significance", 0)) < significance_threshold:
            continue

        # scan all Description fields under this record
        for desc in _find_descriptions(rec):
            # normalize to a single string
            text = ""
            if isinstance(desc, dict):
                text = " ".join(str(v) for v in desc.values())
            else:
                text = str(desc)
            if any(kw in text.lower() for kw in kws):
                # grab the gene name from geneProperties
                name = rec.get("geneProperties", {}).get("gene_name")
                if name:
                    print("-----------------------------------------------------------------")
                    print(extract_en_names(rec))


                    result.extend(extract_en_names(rec) )
                break

    return result


# Example:
if __name__ == "__main__":
    genes = filter_gene_names_by_description("b_sil/mgs.usergenes_Wes_3687.json")
    print("Genes meeting significance ≥ 0.8 with neurological descriptions:", genes)
    print( ', '.join(genes) )
