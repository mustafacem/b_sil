import json
import os
from typing import Any, Dict, List, Set, Union


def find_phenotypes_for_inheritance_genotype(json_file_path: str) -> List[str]:
    """
    Reads JSON data, identifies gene records containing variants matching specific
    inheritance + genotype criteria, and returns unique associated phenotypes.

    Criteria
    --------
    • Inheritance **Autosomal recessive**  → genotype must be **Homozygous variant**  
    • Inheritance **Autosomal dominant**  → genotype may be **Homozygous variant** OR **Heterozygous**

    Parameters
    ----------
    json_file_path : str
        Path to the local JSON file.

    Returns
    -------
    list[str]
        Alphabetically-sorted unique *PhenotypeName* values for all variants that meet
        the rules above.  If no matches (or an error), an empty list is returned.
    """
    matching_phenotype_names: Set[str] = set()

    # ── read the file ───────────────────────────────────────────────────────────
    if not os.path.exists(json_file_path):
        print(f"Error: File not found: {json_file_path}")
        return []
    if not os.access(json_file_path, os.R_OK):
        print(f"Error: Cannot read file (permission denied): {json_file_path}")
        return []

    try:
        with open(json_file_path, "r", encoding="utf-8") as fh:
            data_list: Union[List[Dict[str, Any]], Dict[str, Any]] = json.load(fh)
    except (IOError, OSError) as e:
        print(f"Error reading file {json_file_path}: {e}")
        return []
    except json.JSONDecodeError as e:
        print(f"Error: invalid JSON in {json_file_path} – {e}")
        return []

    # normalise: wrap single-dict JSON into a list so the same loops work
    if isinstance(data_list, dict):
        data_list = [data_list]
    if not isinstance(data_list, list):
        print(f"Error: top-level JSON is neither a list nor a dict.")
        return []

    # ── main walk over records/variants ────────────────────────────────────────
    for gene_record in data_list:
        parent_oid = (gene_record.get("_id") or {}).get("$oid", "Unknown OID")
        variants = gene_record.get("Variants", [])
        if not isinstance(variants, list):
            continue

        for variant in variants:
            if not isinstance(variant, dict):
                continue

            genotype = variant.get("Genotype")
            if not isinstance(genotype, str):
                continue

            variant_id = variant.get("VariantID")        # ← grab ID once here
            pheno_inh_list = variant.get("pheno_inh", [])
            if not isinstance(pheno_inh_list, list):
                continue

            for pheno_info in pheno_inh_list:
                if not isinstance(pheno_info, dict):
                    continue

                inheritance = pheno_info.get("Inheritance")
                phenotype_name = pheno_info.get("PhenotypeName")

                if not (isinstance(inheritance, str)
                        and isinstance(phenotype_name, str)
                        and phenotype_name):
                    continue

                match_found = False

                if inheritance == "Autosomal recessive":
                    match_found = genotype == "Homozygous variant"

                elif inheritance == "Autosomal dominant":
                    match_found = genotype in ("Homozygous variant", "Heterozygous")

                if match_found:
                    reason = f"{inheritance[:2]} + {genotype.split()[0]}"
                    print(f"Match: '{phenotype_name}' ({reason}) – gene OID {parent_oid}, variant {variant_id}")
                    matching_phenotype_names.add(phenotype_name)

    return sorted(matching_phenotype_names)


# ── example usage ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    file_path = "mgs.usergenes_Wes_3687.json"          # adjust to your file
    phenotypes = find_phenotypes_for_inheritance_genotype(file_path)

    print("\nPhenotype names that meet the criteria "
          "(AR+Homozygous or AD+Homozygous/Heterozygous):")
    if phenotypes:
        print(", ".join(phenotypes))
    else:
        print("No matching phenotypes found.")
