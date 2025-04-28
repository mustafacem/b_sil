import json
import os
from typing import Any, Dict, List, Set, Union # Added Set for type hint

def find_phenotypes_for_recessive_homozygous(json_file_path: str) -> List[str]:
  """
  Reads JSON data from a local file path, identifies gene records
  containing at least one variant with 'Autosomal recessive' inheritance
  and 'Homozygous variant' genotype, and returns the unique 'PhenotypeName'
  values associated with those matching variants.

  Args:
    json_file_path: A string containing the path to the local JSON file.

  Returns:
    A list of unique 'PhenotypeName' strings for variants that match
    the criteria (Homozygous, Autosomal Recessive). Returns an empty list
    if no matches are found or if there's an error reading/parsing the file.
  """
  matching_phenotype_names: Set[str] = set() # Use a set to store unique phenotype names
  data_list: Union[List[Dict[str, Any]], Dict[str, Any], None] = None # More specific type hint

  # --- Read data from local file path ---
  if not os.path.exists(json_file_path):
      print(f"Error: File not found at path: {json_file_path}")
      return []
  if not os.access(json_file_path, os.R_OK):
      print(f"Error: Permission denied when trying to read file: {json_file_path}")
      return []

  try:
    with open(json_file_path, 'r', encoding='utf-8') as f:
      try:
          data_list = json.load(f)
      except json.JSONDecodeError as e:
          print(f"Error: Could not decode JSON from file: {json_file_path}")
          print(f"JSONDecodeError: {e}")
          return [] # Return empty list

  except IOError as e:
      print(f"Error reading file {json_file_path}: {e}")
      return []
  except Exception as e:
    print(f"An unexpected error occurred while reading/parsing the file: {e}")
    return []

  # --- Process the loaded data ---
  if not isinstance(data_list, list):
      print(f"Warning: JSON data from {json_file_path} is not a list. Trying to process as single object.")
      if isinstance(data_list, dict):
          data_list = [data_list] # Wrap single dict in a list
      else:
          print(f"Error: Data from {json_file_path} is neither a list nor a dictionary that can be processed.")
          return []

  # --- Iterate through gene records ---
  for gene_record in data_list:
    # Keep track of parent OID for potential debugging/logging, but don't collect it
    parent_oid_str = "Unknown OID" # Default value
    parent_id_dict = gene_record.get("_id")
    if isinstance(parent_id_dict, dict):
        parent_oid = parent_id_dict.get("$oid")
        if isinstance(parent_oid, str):
            parent_oid_str = parent_oid # Store for logging if needed

    # Safely get the Variants list
    variants = gene_record.get("Variants", [])
    if not isinstance(variants, list):
        # print(f"Warning: 'Variants' field is not a list in record with OID {parent_oid_str}. Skipping.")
        continue # Skip if Variants is not a list

    # --- Iterate through variants for this gene record ---
    for variant in variants:
      if not isinstance(variant, dict):
          # print(f"Warning: Found non-dictionary item in 'Variants' list for record {parent_oid_str}. Skipping item.")
          continue # Skip non-dict items in Variants list

      # Safely get genotype
      genotype = variant.get("Genotype")
      if not isinstance(genotype, str) or genotype != "Homozygous variant":
          continue # Skip if not homozygous or genotype missing/invalid

      # Safely get pheno_inh list
      pheno_inh_list = variant.get("pheno_inh", [])
      if not isinstance(pheno_inh_list, list):
          # print(f"Warning: 'pheno_inh' is not a list for a variant in record {parent_oid_str}. Skipping variant.")
          continue # Skip if pheno_inh is not a list

      # --- Check inheritance within pheno_inh ---
      for pheno_info in pheno_inh_list:
        if not isinstance(pheno_info, dict):
            # print(f"Warning: Found non-dictionary item in 'pheno_inh' list for a variant in record {parent_oid_str}. Skipping item.")
            continue # Skip non-dict items in pheno_inh list

        inheritance = pheno_info.get("Inheritance")
        # Check if inheritance is 'Autosomal recessive'
        if isinstance(inheritance, str) and inheritance == "Autosomal recessive":
          # Match found! Get the PhenotypeName
          phenotype_name = pheno_info.get("PhenotypeName")
          if isinstance(phenotype_name, str) and phenotype_name: # Ensure it's a non-empty string
            print(f"Found match: Phenotype='{phenotype_name}' (in Gene Record OID: {parent_oid_str})")
            matching_phenotype_names.add(phenotype_name)
            # No need to break here if a variant can have multiple matching phenotypes
          # else:
            # print(f"Warning: Match found (Homozygous, AR) but 'PhenotypeName' missing or invalid in record {parent_oid_str}.")

        # Removed the inner break: continue checking other pheno_inh entries for this variant

      # Removed the outer break: continue checking other variants in this gene record

  # Convert the set of unique phenotype names to a list before returning
  return sorted(list(matching_phenotype_names)) # Return sorted list for consistent order


# --- Example Usage ---
# Create a dummy directory and file for testing if needed
# Adjust the path as necessary for your actual file


# --- Call the function and print results ---
phenotype_names = find_phenotypes_for_recessive_homozygous("b_sil/mgs.usergenes_Wes_3687.json")

# Check if the function returned a list (it should always, even if empty)
if isinstance(phenotype_names, list):
    print("\nPhenotype Names associated with Homozygous variants with Autosomal Recessive inheritance:")
    if phenotype_names:
        print(( ', '.join(phenotype_names) ))
    else:
        print("No matching phenotypes found.")
# No need for 'is not None' check as the function guarantees a list return on success/handled errors