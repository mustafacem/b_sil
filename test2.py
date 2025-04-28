import json, os
from typing import Any, Dict, List, Set, Union


def _has_pathogenic_label(variant: Dict[str, Any]) -> bool:
    """
    True if *any* Final_label for this variant (top-level or in pheno_inh)
    contains the word “pathogenic” (case-insensitive).
    """
    # ── top-level Final_label ──────────────────────────────────────────────
    top = variant.get("Final_label", [])
    if isinstance(top, str):
        top = [top]
    if any("pathogenic" in str(lbl).lower() for lbl in top):
        return True

    # ── per-phenotype Final_label inside pheno_inh ─────────────────────────
    for item in variant.get("pheno_inh", []):
        lbl = item.get("Final_label")
        if isinstance(lbl, str) and "pathogenic" in lbl.lower():
            return True
        if isinstance(lbl, list) and any("pathogenic" in str(x).lower() for x in lbl):
            return True

    return False


def find_phenotypes_for_inheritance_genotype(json_file_path: str) -> List[str]:
    """
    Return VariantIDs that satisfy BOTH conditions:

      1. Inheritance/Genotype rule  
           • Autosomal recessive → Homozygous variant  
           • Autosomal dominant → Homozygous variant *or* Heterozygous
      2. At least one Final_label contains the word “Pathogenic”
         (e.g. “Likely pathogenic”, “Pathogenic*”, …)

    On error or if nothing matches, an empty list is returned.
    """
    matching_variant_ids: Set[str] = set()

    # ── read JSON file ─────────────────────────────────────────────────────
    if not (os.path.exists(json_file_path) and os.access(json_file_path, os.R_OK)):
        print(f"Error: cannot open {json_file_path}")
        return []

    try:
        with open(json_file_path, "r", encoding="utf-8") as fh:
            data: Union[List[Dict[str, Any]], Dict[str, Any]] = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        print(f"Error reading {json_file_path}: {e}")
        return []

    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        print("Error: top-level JSON must be a list or dict.")
        return []

    # ── walk the file ──────────────────────────────────────────────────────
    for gene_record in data:
        for variant in gene_record.get("Variants", []):
            if not isinstance(variant, dict):
                continue

            genotype = variant.get("Genotype", "")
            if not isinstance(genotype, str):
                continue

            # check every phenotype-inheritance pair stored in pheno_inh
            for ph_info in variant.get("pheno_inh", []):
                inh  = ph_info.get("Inheritance")
                if not isinstance(inh, str):
                    continue

                rule_ok = (
                    (inh == "Autosomal recessive" and genotype == "Homozygous variant")
                    or
                    (inh == "Autosomal dominant"  and genotype in ("Homozygous variant", "Heterozygous"))
                )

                if rule_ok and _has_pathogenic_label(variant):
                    matching_variant_ids.add(variant.get("VariantID", ""))

    return sorted(v for v in matching_variant_ids if v)   # drop empty strings


# ── demo ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    path = "mgs.usergenes_Wes_3687.json"
    variant_ids = find_phenotypes_for_inheritance_genotype(path)

    print("\nVariantIDs that meet inheritance/genotype rules *and* have a "
          "'Pathogenic' Final_label:")
    print(", ".join(variant_ids) if variant_ids else "None found.")
