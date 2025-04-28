import streamlit as st
from openai import OpenAI
import pandas as pd
import json
from datetime import datetime
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
import os 
from typing import Any, Dict, List, Optional, Union
from pathlib import Path

from dotenv import load_dotenv
load_dotenv()



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






def generate_response_com(prompt: str) -> str:
    """
    Send a single-user prompt to the Chat API and return the assistant's reply.
    """
    client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key=os.getenv("OPENR", "")
    )
    # Wrap the prompt in the expected messages format
    messages = [
        {"role": "user", "content": prompt}
    ]
    completion = client.chat.completions.create(
        extra_headers={
            "HTTP-Referer": "<YOUR_SITE_URL>",
            "X-Title": "<YOUR_SITE_NAME>",
        },
        model="google/gemini-2.5-pro-preview-03-25",
        messages=messages
    )
    # Extract and return the assistant’s reply
    return completion.choices[0].message.content.strip()

def check_clinvar_label(json_path: str, variant_id: str) -> str:
    """
    Checks the ClinVar label for a given variant.

    Args:
        json_path: Path to the JSON file.
        variant_id: The VariantID to look up.

    Returns:
        "yes" if the ClinVar annotation is present and not "None",
        "no" if ClinVar == "None" or missing,
        raises KeyError if the variant_id isn't found.
    """
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    for var in data.get("variants", []):
        if var.get("VariantID") == variant_id:
            clinvar = (
                var
                .get("annotations", {})
                .get("Variant databases", {})
                .get("ClinVar")
            )
            # If clinvar is explicitly "None" or not set → "no"
            if clinvar is None or clinvar == "None":
                return "no"
            # Otherwise → "yes"
            return "yes"

    # Variant not found
    raise KeyError(f"VariantID '{variant_id}' not found in {json_path}")

def get_final_score(variant_id: str, json_path: str) -> float:
    """
    Load JSON from `json_path`, search all variants for `variant_id`,
    and return its Final_score as a float.
    
    Raises:
        ValueError: if the JSON structure is unexpected,
                    if the variant isn’t found,
                    or if Final_score is missing/invalid.
    """
    # 1. Load the JSON file
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # 2. Gather all variants into one list
    variants: List[Dict[str, Any]] = []
    def extract(rec: Dict[str, Any]):
        for key in ('variants', 'Variants'):
            lst = rec.get(key)
            if isinstance(lst, list):
                variants.extend(lst)

    if isinstance(data, dict):
        extract(data)
    elif isinstance(data, list):
        for record in data:
            if isinstance(record, dict):
                extract(record)
    else:
        raise ValueError("JSON must be a dict or a list of dicts")

    # 3. Find the matching variant and return its Final_score
    for var in variants:
        if var.get('VariantID') == variant_id:
            score = var.get('Final_score')
            if score is None:
                raise ValueError(f"'Final_score' missing for variant {variant_id}")
            try:
                return float(score)
            except (TypeError, ValueError):
                raise ValueError(f"Invalid Final_score value: {score!r}")

    # 4. Not found
    raise ValueError(f"VariantID '{variant_id}' not found")


def find_unique_phenotypes_by_gene(json_path: str, gene_id: str) -> List[Dict[str, Any]]:
    """
    Load JSON (single dict or list of dicts) and return unique phenotype entries
    where the given gene_id appears in the Genes section.

    Args:
        json_path: Path to the JSON file.
        gene_id: The gene symbol to search for (e.g. "VAMP7").

    Returns:
        A list of unique phenotype dicts matching that gene.
    """
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Normalize to a list of records
    records = data if isinstance(data, list) else [data]

    seen_ids = set()
    unique_phenos: List[Dict[str, Any]] = []

    for record in records:
        for variant in record.get("Variants", []):
            for pheno in variant.get("phenotypes", []):
                # Determine a unique key for this phenotype
                # Prefer the phenotypeID field if present; otherwise fall back to internal $oid
                pheno_id = (
                    pheno.get("phenotypeID")
                    or pheno.get("_id", {}).get("$oid")
                    or None
                )
                if not pheno_id or pheno_id in seen_ids:
                    continue

                # Check if gene_id is in any Genes mapping
                for gene_map in pheno.get("Genes", []):
                    if gene_id in gene_map:
                        seen_ids.add(pheno_id)
                        unique_phenos.append(pheno)
                        break  # stop checking this pheno once matched

    return unique_phenos

def get_carrier_message(variant_id: str, json_path: str) -> str:
    """
    Load JSON from `json_path`, search all variants for `variant_id`,
    and return:
      - "Your children will be carriers"   if Genotype == "Homozygous variant"
      - "Your children might be carriers" if Genotype == "Heterozygous"
      - an informative message otherwise.
    """
    # 1. Load the JSON
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # 2. Gather all variant dicts into a single list
    variants: List[Dict[str, Any]] = []
    
    def extract_from_record(rec: Dict[str, Any]):
        for key in ('variants', 'Variants'):
            vals = rec.get(key)
            if isinstance(vals, list):
                variants.extend(vals)

    if isinstance(data, dict):
        extract_from_record(data)
    elif isinstance(data, list):
        for record in data:
            if isinstance(record, dict):
                extract_from_record(record)
    else:
        raise ValueError("JSON must be a dict or a list of dicts")

    # 3. Find the matching VariantID
    for var in variants:
        if var.get('VariantID') == variant_id:
            geno = var.get('Genotype', '')
            if geno == 'Homozygous variant':
                return 'Your children will be carriers'
            if geno == 'Heterozygous':
                return 'Your children might be carriers'
            return f"Unrecognized genotype: '{geno}'"

    # 4. Not found
    return f"VariantID '{variant_id}' not found"  


def get_phenotype_details(
    json_path: str,
    oids: Union[str, List[str]]
) -> Dict[str, Dict[str, Any]]:
    """
    Load JSON from `json_path`, recursively scan for phenotype objects, and
    return a dict mapping each requested oid to its details:
      - Name (English)
      - Inheritance (English list)
      - Symptoms (English list)
      - Description (English)
      - AgeOfOnset (English list)

    Args:
        json_path: Path to your JSON file (which may be a dict or a list at top level).
        oids: A single oid string or a list of oid strings to look up.

    Returns:
        A dict where each key is one of the requested oids (if found), and each
        value is another dict with keys "Name", "Inheritance", "Symptoms",
        "Description", and "AgeOfOnset".
    """
    # Load the JSON file
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)

    # Normalize oids to a list
    if isinstance(oids, str):
        oids = [oids]

    # This will collect details for any matching oid
    details: Dict[str, Dict[str, Any]] = {}

    def recurse(obj: Any):
        if isinstance(obj, dict):
            oid = obj.get("_id", {}).get("$oid")
            # If this dict has an oid we care about, extract its fields
            if oid in oids:
                details[oid] = {
                    "Name": obj.get("Name", {}).get("en"),
                    "Inheritance": obj.get("Inheritance", {}).get("en", []),
                    "Symptoms": obj.get("Symptoms", {}).get("en", []),
                    "Description": obj.get("Description", {}).get("en"),
                    "AgeOfOnset": obj.get("AgeOfOnset", {}).get("en", [])
                }
            # Recurse into all values
            for v in obj.values():
                recurse(v)
        elif isinstance(obj, list):
            for item in obj:
                recurse(item)

    # Build the map
    recurse(data)

    return details


def find_significant_phenotypes_from_file(json_file_path, threshold=0.7):
    """
    Reads JSON data from a file, parses it, and returns phenotype OIDs
    for entries with significance >= threshold.

    Args:
        json_file_path: Path to the JSON file containing a single JSON object
                        or a list of JSON objects.
        threshold: The minimum significance value (inclusive) to filter by.
                   Defaults to 0.7.

    Returns:
        A list of phenotype '$oid' strings that meet the significance criteria.
        Returns an empty list if the file cannot be read, is not valid JSON,
        no matches are found, or the JSON structure is unexpected.
    """
    significant_oids = []
    data = None

    # 1. Read and parse the JSON file
    try:
        with open(json_file_path, 'r', encoding='utf-8') as f: # Use utf-8 encoding
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File not found at path: {json_file_path}")
        return []
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from file {json_file_path}: {e}")
        return []
    except IOError as e: # Catch other potential file reading errors
        print(f"Error reading file {json_file_path}: {e}")
        return []
    except Exception as e: # Catch unexpected errors during file handling/parsing
        print(f"An unexpected error occurred while processing {json_file_path}: {e}")
        return []


    # 2. Determine if we have a single object or a list from the parsed data
    items_to_check = []
    if isinstance(data, dict):
        # Wrap single dictionary in a list to use the same loop logic
        items_to_check = [data]
    elif isinstance(data, list):
        items_to_check = data
    else:
        # This case might occur if the JSON file contained something other than an object or array
        print(f"Warning: Parsed data from {json_file_path} is not a dictionary or list: {type(data)}")
        return []

    # 3. Iterate and check significance (same logic as before)
    for item in items_to_check:
        # Ensure the item is actually a dictionary before proceeding
        if not isinstance(item, dict):
            print(f"Warning: Skipping non-dictionary item found in list from file {json_file_path}: {item}")
            continue

        # Safely get significance and phenotype info using .get()
        significance = item.get('significance')
        phenotype_info = item.get('phenotype')

        # Check if significance exists, is a number, and meets the threshold
        if significance is not None and isinstance(significance, (int, float)) and significance >= threshold:
            # Check if phenotype info exists and is a dictionary
            if phenotype_info and isinstance(phenotype_info, dict):
                phenotype_oid = phenotype_info.get('$oid')
                # Add the oid if it exists
                if phenotype_oid and isinstance(phenotype_oid, str): # Ensure oid is a string
                    significant_oids.append(phenotype_oid)
                # Optionally add a warning if $oid is missing or not a string
                # else:
                #    print(f"Warning: Found significant item with phenotype dict but missing or invalid '$oid': {item}")
            # Optionally add a warning if phenotype is missing/invalid for a significant item
            # else:
            #    print(f"Warning: Found significant item but missing or invalid 'phenotype' dictionary: {item}")

    return significant_oids


def get_names_by_oid(json_path: str,
                     oids: Union[str, List[str]]) -> List[str]:
    """
    Load JSON from `json_path`, recursively scan for any objects
    with {"_id": {"$oid": ...}, "Name": {"en": ...}}, and return the
    English names for the given oid(s).

    Args:
        json_path: Path to the JSON file (which may contain a dict or a list).
        oids: A single oid string or a list of oid strings to look up.

    Returns:
        A list of Name["en"] values for each oid found (in the order requested).
    """
    # Load the JSON file
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)

    # Normalize oids to a list
    if isinstance(oids, str):
        oids = [oids]

    # Collector for oid → name_en
    oid_to_name: Dict[str, str] = {}

    def recurse(obj: Any):
        if isinstance(obj, dict):
            oid = obj.get("_id", {}).get("$oid")
            name_en = obj.get("Name", {}).get("en")
            if oid and name_en:
                oid_to_name[oid] = name_en
            # Recurse into all values
            for v in obj.values():
                recurse(v)
        elif isinstance(obj, list):
            for item in obj:
                recurse(item)

    # Build the map
    recurse(data)

    # Return names for the requested oids (skip any not found)
    return [oid_to_name[oid] for oid in oids if oid in oid_to_name]
                         
def generate_pdf(
    json_file: str | Path = "mgs.userphenotypes_Wes_3687.json",
    pdf_filename: str | Path = "genetic_analysis_report.pdf",
    patient_name: str = "",
    patient_gender: str = "Male",
    final_score_threshold: float = 0.9,
    wide_report: bool = False,
) -> None:
    """
    Build a PDF listing only those variants whose Final_score > threshold.
    * Works whether the JSON root is a dict or list.
    * Ensures each VariantID appears only once (first occurrence kept).
    * Wide layout now **omits** the former “Extra1/Extra2” columns.
    """
    # ── 1. load & normalise ───────────────────────────────────────────────
    with open(json_file, encoding="utf-8") as fh:
        raw = json.load(fh)

    records: list[dict] = raw if isinstance(raw, list) else [raw]

    variants_by_id: dict[str, dict] = {}       # keeps insertion order
    information_en = "No additional information provided."

    for rec in records:
        if not isinstance(rec, dict):
            continue

        for v in rec.get("variants", []):      # sample file uses lowercase key
            vid = v.get("VariantID")
            if vid and vid not in variants_by_id:
                variants_by_id[vid] = v

        if "information" in rec and "en" in rec["information"]:
            information_en = rec["information"]["en"]

    variants = list(variants_by_id.values())

    # ── 2. filter by Final_score ──────────────────────────────────────────
    def score_is_high(v: dict) -> bool:
        val = v.get("Final_score")
        try:
            return float(val) > final_score_threshold
        except (TypeError, ValueError):
            return False

    filtered = [v for v in variants if score_is_high(v)]

    # ── 3. build the PDF ──────────────────────────────────────────────────
    current_date = datetime.now().strftime("%d.%m.%Y")
    doc    = SimpleDocTemplate(pdf_filename, pagesize=letter)
    styles = getSampleStyleSheet()
    story  = []

    # Header
    story.append(Paragraph("<b>MGS-AI Chat Report Sample</b>", styles["Title"]))
    story.append(Spacer(1, 12))
    story.append(Paragraph("<b>Genetic Analysis Report</b>", styles["Heading1"]))
    story.append(Spacer(1, 24))

    # Patient information
    story.append(Paragraph("<b>Patient Information</b>", styles["Heading2"]))
    story.append(Spacer(1, 12))
    if patient_name.strip():
        story.append(Paragraph(f"Name: {patient_name.strip()}", styles["Normal"]))
        story.append(Spacer(1, 6))
    story.append(Paragraph(f"Gender: {patient_gender}", styles["Normal"]))
    story.append(Spacer(1, 6))
    story.append(Paragraph(f"Date: {current_date}", styles["Normal"]))
    story.append(Spacer(1, 24))

    # Summary
    story.append(Paragraph("<b>Summary</b>", styles["Heading2"]))
    story.append(Spacer(1, 12))
    summary_html = (
        "This report summarizes the results of automated variant interpretation for the patient. "
        f"Only variants with a Final Score greater than <b>{final_score_threshold:.2f}</b> are displayed. "
        "Variants have been analyzed based on phenotype matching, variant pathogenicity, inheritance pattern, "
        "and supporting literature.<br/><br/>"
        f"<i>Additional Info:</i> {information_en}"
    )
    story.append(Paragraph(summary_html, styles["Normal"]))
    story.append(Spacer(1, 24))

    # Variants table
    story.append(Paragraph("<b>Variants Data</b>", styles["Heading2"]))
    story.append(Spacer(1, 12))

    header = (
        ["Variant", "Type", "Genotype", "Gene", "Label"]         # ← wide header without extras
        if wide_report
        else ["Variant", "Type", "Genotype", "Gene", "Phenotype", "Inheritance Type", "Classification"]
    )
    table_data = [header]

    for v in filtered:
        if wide_report:
            row = [
                v.get("VariantID", "N/A"),
                v.get("Type", "N/A"),
                v.get("Genotype", "N/A"),
                v.get("Gene", "N/A"),
                v.get("Label", "N/A"),
            ]
        else:
            row = [
                v.get("VariantID", "N/A"),
                v.get("Type", "N/A"),
                v.get("Genotype", "N/A"),
                v.get("Gene", "N/A"),
                "*",  # placeholder until phenotype data is available
                "*",  # placeholder until inheritance data is available
                v.get("Label", "N/A"),
            ]
        table_data.append(row)

    vt = Table(table_data, hAlign="CENTER")
    vt.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 1, colors.black),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
                ("TOPPADDING", (0, 0), (-1, 0), 8),
            ]
        )
    )
    story.append(vt)
    story.append(Spacer(1, 24))

    # Variant details
    story.append(Paragraph("<b>Variant Details</b>", styles["Heading2"]))
    story.append(Spacer(1, 12))
    for v in filtered:
        vid = v.get("VariantID", "N/A")
        story.append(Paragraph(vid, styles["Heading3"]))
        story.append(Spacer(1, 6))

        ann   = v.get("annotations", {})
        db    = ann.get("Variant databases", {})
        freq  = ann.get("Population allele frequency", {})
        pred  = ann.get("Pathogenicity predictions", {})

        details = [
            ["ClinVar", "GnomAD", "PHACTboost", "AlphaMissense"],
            [
                db.get("ClinVar", "N/A"),
                freq.get("GnomAD joint allele frequency", "N/A"),
                pred.get("PHACTboost", "N/A"),
                pred.get("AlphaMissense", "N/A"),
            ],
        ]
        dt = Table(details, hAlign="CENTER")
        dt.setStyle(TableStyle([("ALIGN", (0, 0), (-1, -1), "CENTER")]))
        story.append(dt)
        story.append(Spacer(1, 12))

    # Results & notes
    story.append(Paragraph("<b>Results</b>", styles["Heading2"]))
    story.append(Spacer(1, 12))
    story.append(
        Paragraph(
            "The genetic variants listed above were identified and characterized based on the provided data. "
            "Further clinical correlation is advised given the variant classifications.",
            styles["Normal"],
        )
    )
    story.append(Spacer(1, 24))

    story.append(Paragraph("<b>Notes</b>", styles["Heading2"]))
    story.append(Spacer(1, 12))
    story.append(
        Paragraph(
            "Reviewed by Dr. A. Smith. Recommend further clinical follow-up for significant findings.",
            styles["Normal"],
        )
    )

    doc.build(story)
    print(
        f"PDF '{pdf_filename}' generated with {len(filtered)} variant(s) "
        f"having Final Score > {final_score_threshold}.",
    )



# ------------------------------------------------------------------------------

df_questions = pd.read_excel("AI chat - potential questions.xlsx")

def generate_response(messages):
    client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key=os.getenv("OPENR", "")
    )
    completion = client.chat.completions.create(
        extra_headers={
            "HTTP-Referer": "<YOUR_SITE_URL>",
            "X-Title": "<YOUR_SITE_NAME>",
        },
        model="google/gemini-2.5-pro-preview-03-25",
        messages=messages
    )
    return completion.choices[0].message.content

def match_user_question(user_query, df, llm_client):
    """
    Returns the matching row (as a dict) from df based on the user query,
    or returns None if no match is found.
    """
    possible_questions = [q for q in df["Questions"]]
    system_prompt = (
        "You are a helper AI. The user asks a question, and you have a list of known questions or tasks. "
        "Pick which one from the list below best matches the user's question. "
        "If none match, return 'NONE'.\n\n"
        "List of known questions:\n"
    )
    for i, q in enumerate(possible_questions):
        system_prompt += f"{i+1}. {q}\n"
    system_prompt += "\nRespond ONLY with the number of the best matching question or 'NONE'."

    mini_messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_query}
    ]
    completion = llm_client.chat.completions.create(
        extra_headers={"HTTP-Referer": "<YOUR_SITE_URL>", "X-Title": "<YOUR_SITE_NAME>"},
        model="google/gemini-2.5-pro-preview-03-25",
        messages=mini_messages
    )
    match_response = completion.choices[0].message.content.strip().upper()
    if "NONE" in match_response:
        return None
    try:
        matched_index = int("".join([c for c in match_response if c.isdigit()])) - 1
        if matched_index < 0 or matched_index >= len(df):
            return None
        return df.iloc[matched_index].to_dict()
    except:
        return None

def parse_report_details(details_text):
    """
    Infers the patient's name and the report type (Normal or Wide) from text.
    Returns (patient_name, report_type).
    """
    prompt = (
        "You are an expert assistant specialized in processing natural language report requests. "
        "Extract the patient's name and the report type from the following text. The report type must be either "
        "'Normal' or 'Wide'. If the patient's name is not specified, return an empty string for the patient name. "
        "Return exactly valid JSON with keys: 'patient_name' and 'report_type'."
        "\n\nInput: " + details_text
    )
    messages = [{"role": "system", "content": prompt}]
    llm_response = generate_response(messages)
    try:
        result = json.loads(llm_response)
        patient_name = result.get("patient_name", "").strip()
        report_type = result.get("report_type", "Normal").strip().capitalize()
        if report_type not in ("Normal", "Wide"):
            report_type = "Normal"
        return (patient_name, report_type)
    except Exception:
        return ("", "Normal")


def get_phenotypes_for_risk():
    return ["Phenotype A", "Phenotype B", "Phenotype C"]


def main():
    # Page configuration
    st.set_page_config(page_title="MGS Personal Chat", layout="centered")

    # Optional styling
    st.markdown(
        """
        <style>
        /* Center container */
        .main > div {
            max-width: 700px !important;
            margin: 0 auto;
        }
        .stChatMessage {
            margin-bottom: 1rem;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

    st.title("MGS Personal Chat")

    # Initialize session state variables
    if "messages" not in st.session_state:
        st.session_state["messages"] = []
    if "await_report_details" not in st.session_state:
        st.session_state["await_report_details"] = False
    if "report_generated" not in st.session_state:
        st.session_state["report_generated"] = False

    # Capture user input
    user_input = st.chat_input("Type your message and press Enter...")

    if user_input:
        # If awaiting PDF report details
        if st.session_state["await_report_details"]:
            patient_name, report_type = parse_report_details(user_input)
            if report_type not in ("Normal", "Wide"):
                st.session_state["messages"].append({
                    "role": "assistant",
                    "content": (
                        "I couldn't determine the report type. "
                        "Please specify 'Normal' or 'Wide', e.g. 'John Doe | Normal'."
                    )
                })
            else:
                # Generate PDF based on report type
                if report_type == "Wide":
                    generate_pdf(
                        pdf_filename="genetic_analysis_report.pdf",
                        patient_name=patient_name,
                        patient_gender="Male",
                        final_score_threshold=0.9,
                        wide_report=True
                    )
                else:
                    generate_pdf(
                        pdf_filename="genetic_analysis_report.pdf",
                        patient_name=patient_name,
                        patient_gender="Male",
                        final_score_threshold=0.45,
                        wide_report=False
                    )

                # Notify user and set download flag
                st.session_state["messages"].append({
                    "role": "assistant",
                    "content": f"✅ PDF generated as a {report_type.lower()} report."
                })
                st.session_state["report_generated"] = True
                st.session_state["await_report_details"] = False

        else:
            # Normal chat scenario
            st.session_state["messages"].append({"role": "user", "content": user_input})

            # Attempt to match known questions
            match_client = OpenAI(
                base_url="https://openrouter.ai/api/v1",
                api_key=os.getenv("OPENR", "")
            )
            matched_row = match_user_question(user_input, df_questions, match_client)

            if matched_row:
                matched_answer = matched_row.get("Answers", "")
                user_question = matched_answer.lower()
                row_question  = matched_row.get("Questions", "").lower()

                if ("create a report" in user_question
                    or row_question == "create me a report"):
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (
                            "Let's create your report. Please provide the details in one of the following ways:\n\n"
                            "1. In the format: 'Patient Name | Normal' or 'Patient Name | Wide'\n"
                            "2. Or simply type your request naturally (e.g., 'generate a wide report for John Doe')."
                        )
                    })
                    st.session_state["await_report_details"] = True

                elif ("what are the phenotypes i have risk for?" in user_question
                    or row_question == "what are the phenotypes i have risk for?"):
                    phens_to_return = find_significant_phenotypes_from_file(
                        "mgs.userphenotypes_Wes_3687.json", 0.7
                    )
                    print("phens_to_return", phens_to_return)
                    phens_to_print = get_names_by_oid(
                        "mgs.usergenes_Wes_3687.json",
                        phens_to_return
                    )
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": f"you  have risk for these phenotypes :{phens_to_print}" 
                    })

                elif ("what are the symptoms of this phenotype?" in user_question
                    or row_question == "what are the symptoms of this phenotype?"):
                      #beni sonra değişlitir  dyncamic olmalı 

                    """
                    phens_to_return = find_significant_phenotypes_from_file(
                        "pheno_input.json", 0.7
                    )
                    print("phens_to_return", phens_to_return)
                    details = get_phenotype_details(
                        "mgs.usergenes_Wes_3687.json",
                        phens_to_return
                    )

                    """
                    output_parts = []
                    
                    data_ofmsg = st.session_state["messages"]
                    print("details", data_ofmsg)

                    assistant_content = next(
                        (message['content'] for message in reversed(data_ofmsg) if message.get('role') == 'assistant'),
                        None  # Default value if no assistant message is found
                    )
                    print("onceki", assistant_content)                   
                    phens_to_return = generate_response_com( "find a return phenotypeID which looks like : 674f14f21bb9c321872dc66d , from this text return nothing else:" + assistant_content)
                    print("9999999",  phens_to_return)
                    phens_to_return = [phens_to_return]
                    details = get_phenotype_details(
                        "mgs.usergenes_Wes_3687.json",
                        phens_to_return
                    )
                   




                    for oid, info in details.items():
                        name = info.get("Name", "<no name>")
                        symptoms = info.get("Symptoms", [])
                        symptoms_str = ", ".join(symptoms) or "None"
                        # one line per phenotype
                        output_parts.append(f"{name} -> {symptoms_str}")

                    output = "\n".join(output_parts)
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": output
                    })

                elif ("what are my causative/high-impact variants?" in user_question # inhretance yok 680c217dc8cb5cb4c6976cea
                    or row_question == "what are my causative/high-impact  variants?"):
                    print("Matched row:", matched_row)
                    phenotype_names = find_phenotypes_for_inheritance_genotype("mgs.usergenes_Wes_3687.json")
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (f"Your  causative/high-impact  variants are :  {( ', '.join(phenotype_names) )}"  
                        )
                    })

                elif ("which phenotypes related to heart and cardiovascular health are present in my genetic data?" in user_question
                    or row_question == "which phenotypes related to heart and cardiovascular health are present in my genetic data?"): # nasil bilcem neuroloogival riskleri , sentomlarda x kelimnesi olması yeterli mi
                    print("Matched row:", matched_row)
                    heart_gens = filter_gene_names_by_description("mgs.usergenes_Wes_3687.json",0.7,["hearth", "cardiovascular","heart","cardiac"])
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (f"You have genetic predisposition for these phenotypes:  {( ', '.join(heart_gens) )}"  )
                    })
                elif ("what are my genetic risks for neurological disorders?" in user_question # nasil bilcem neuroloogival riskleri , sentomlarda x kelimnesi olması yeterli mi // descrpition neurological neuron brain 
                    or row_question == "what are my genetic risks for neurological disorders?"):
                    print("Matched row:", matched_row)
                    neuron_gens = filter_gene_names_by_description("mgs.usergenes_Wes_3687.json")
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (f"You have genetic predisposition for these phenotypes:  {( ', '.join(neuron_gens) )}"  )
                    })

                elif ("is this variant present in clinvar/gnomad?" in user_question
                    or row_question == "is this variant present in clinvar/gnomad?"):
                    print("Matched row:", matched_row)
                    path = "personalchat/pheno_input.json"     
                    vid  = "X-154835925-C-T"
    
                    yesminomu= check_clinvar_label(path, vid)
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (
                            f"This variant is {'present' if yesminomu else 'not present'} in ClinVar/GnomAD."
                        )
                    })

                elif ("what are the phenotypes i am carrier for?" in user_question
                    or row_question == "what are the phenotypes i am carrier for?"):
                    print("Matched row:", matched_row) #neden ikisi de var, genlerde phenotype inhretance yok 680c217dc8cb5cb4c6976cea , variant içerisinde  otomal recessive ise phenotipe varsa 1 yeter  ikiside sağlanack
                    phenotype_names = find_phenotypes_for_recessive_homozygous("mgs.usergenes_Wes_3687.json")
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (    f"You are carrier for:  {( ', '.join(phenotype_names) )}"             )
                    })

                elif ("what are the phenotypes that this gene is related?" in user_question
                    or row_question == "what are the phenotypes that this gene is related?"):
                    print("Matched row:", matched_row)
                    gene = "VAMP7"

                    genes_that_are_related =  find_unique_phenotypes_by_gene(gene, "mgs.usergenes_Wes_3687.json")
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (f"unique phenotypes related to inputed gene: {genes_that_are_related}"
                        )
                    })

                elif ("why is this variant labelled as pathogenic?" in user_question
                    or row_question == "why is this variant labelled as pathogenic?"): # ClinVar , PHACTboost, AlphaMissense, GnomAD joint allele frequency #phenotypeID
                    output_parts = []

                    print("asdsad row:", user_input)
                    phens_to_return = generate_response_com( "find a return variantID which looks like  random id for example : X-156001632-T-A, from this text return nothing else:" + user_input)
                    print("9999999",  phens_to_return)

                    details = get_variant_annotations(
                        phens_to_return,
                        "mgs.usergenes_Wes_3687.json"
                        
                    )
                    items_str = [f"{key}: {value}" for key, value in details.items()]

                    # Join the list elements with ", "
                    output_str = ", ".join(items_str)

                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": f" Variant {phens_to_return}  has these labels and scores  in databases and pathogenicity predictors:  {output_str}"
                    })

                elif ("is this variant rare or common?" in user_question
                    or row_question == "is this variant rare or common?"):
                    print("Matched row:", matched_row) #asil soru 
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (
                            "Based on allele frequency databases, this variant is considered …"
                        )
                    })

                elif ("what genes contain pathogenic variants in my data?" in user_question
                    or row_question == "what genes contain pathogenic variants in my data?"):
                    print("Matched row:", matched_row)#bulamadım anlamadım
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (
                            "These genes have one or more pathogenic variants in your dataset: …"
                        )
                    })

                elif ("does this variant affect my children’s health?" in user_question
                    or row_question == "does this variant affect my children’s health?"):
                    print("Matched row:", matched_row) # id nerden gelicek?
                    child_info =  (get_carrier_message('X-156001632-T-A', "mgs.usergenes_Wes_3687.json"))

                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (child_info)                        
                    })

                elif ("what are my autosomal recessive risk factors?" in user_question # inhretance yok 680c217dc8cb5cb4c6976cea
                    or row_question == "what are my autosomal recessive risk factors?"):
                    print("Matched row:", matched_row)
                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": (
                            "Your autosomal‐recessive risks include these phenotypes: …"
                        )
                    })

                elif ("how is this variant inherited?" in user_question # realted derken parent bakcan 
                    or row_question == "how is this variant inherited?"):

                    output_parts = []

                    print("asdsad row:", user_input)
                    phens_to_return = generate_response_com( "find a return variantID which looks like  random id for example : X-156001632-T-A, from this text return nothing else:" + user_input)
                    print("9999999",  phens_to_return)


                    path = "mgs.usergenes_Wes_3687.json"          # the file you showed
                    vid  = phens_to_return      
                    dict1 = get_inheritance_for_variant(path, vid)
                    result_dict1 = {key: value[0] for key, value in dict1.items()}              # pick any VariantID present
                    items_str = [f"{value}: {key}" for key, value in result_dict1.items()]

                    output_str = ", ".join(items_str)

                    st.session_state["messages"].append({
                        "role": "assistant",
                        "content": f" Variant {phens_to_return} has inheritence of  {output_str.replace(":", " for ")}"
                    })


                    

                else:
                    # Fallback to returning the matched answer
                    print("mo_nmatch")
                    print("Matched row:", matched_row)
                    print("/////",matched_answer, matched_row.get("Questions", ""))

                    st.session_state["messages"].append({
                        "role": "assistant", 
                        "content": matched_answer
                    })

            else:
                # Fallback to LLM chat
                print("mo_nmatch232")

                llm_reply = generate_response(st.session_state["messages"])
                st.session_state["messages"].append({"role": "assistant", "content": llm_reply})

    # Display the chat history
    for msg in st.session_state["messages"]:
        if msg["role"] == "user":
            with st.chat_message("user"):
                st.markdown(msg["content"])
        else:
            with st.chat_message("assistant"):
                st.markdown(msg["content"])

    # If a report was just generated, show download button
    if st.session_state["report_generated"]:
        with st.chat_message("assistant"):
            with open("genetic_analysis_report.pdf", "rb") as pdf_file:
                st.download_button(
                    label="📥 Download your PDF report",
                    data=pdf_file,
                    file_name="genetic_analysis_report.pdf",
                    mime="application/pdf"
                )
        # Reset the flag so the button shows only once
        st.session_state["report_generated"] = False

if __name__ == "__main__":
    main()
