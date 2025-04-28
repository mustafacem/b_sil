import json
import os
from typing import Any, Dict, List, Set, Union # Added Set for type hint

import json, os
from typing import Any, Dict, List, Set, Union

def _has_pathogenic_label(variant: Dict[str, Any]) -> bool:
    """
    True if *any* Final_label for this variant (top-level or in pheno_inh)
    contains the word “pathogenic” (case-insensitive).
    """
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


def find_phenotypes_for_recessive_homozygous(json_file_path: str) -> List[str]:
    """
    Return unique PhenotypeName strings for variants that are:
    • Genotype == 'Homozygous variant'
    • Have *any* pheno_inh entry with Inheritance == 'Autosomal recessive'
    • Contain the word ‘pathogenic’ in any Final_label
    """
    if not (os.path.exists(json_file_path) and os.access(json_file_path, os.R_OK)):
        return []

    try:
        with open(json_file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return []

    # normalise root into a list
    if isinstance(data, dict):
        data = [data]
    elif not isinstance(data, list):
        return []

    matches: Set[str] = set()

    for rec in data:
        variants = rec.get("Variants", [])
        if not isinstance(variants, list):
            continue

        for var in variants:
            if (
                isinstance(var, dict)
                and var.get("Genotype") == "Homozygous variant"
                and _has_pathogenic_label(var)           # NEW FILTER
            ):
                for ph in var.get("pheno_inh", []):
                    if (
                        isinstance(ph, dict)
                        and ph.get("Inheritance") == "Autosomal recessive"
                    ):
                        name = ph.get("PhenotypeName")
                        if isinstance(name, str) and name:
                            matches.add(name)

    return sorted(matches)

# --- Example Usage ---
# Create a dummy directory and file for testing if needed
# Adjust the path as necessary for your actual file


# --- Call the function and print results ---
phenotype_names = find_phenotypes_for_recessive_homozygous("mgs.usergenes_Wes_3687.json")

# Check if the function returned a list (it should always, even if empty)
if isinstance(phenotype_names, list):
    print("\nPhenotype Names associated with Homozygous variants with Autosomal Recessive inheritance:")
    if phenotype_names:
        print(( ', '.join(phenotype_names) ))
    else:
        print("No matching phenotypes found.")
