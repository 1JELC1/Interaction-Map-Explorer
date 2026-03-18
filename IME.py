"""
Automates CP analysis with Multiwfn and interface extraction between molecular fragments.
"""

import os
import glob
import itertools
import shutil
import multiprocessing
import math
import IME_modules as ime
from pathlib import Path

# Local Multiwfn execution functions

def _detect_num_cpus(ratio: float = 0.7) -> int:
    override = os.getenv("MULTIWFN_THREADS") or os.getenv("MWFN_THREADS")
    if override:
        try:
            n = int(override)
            return max(1, n)
        except Exception:
            pass
    try:
        available = len(os.sched_getaffinity(0))
    except Exception:
        try:
            available = multiprocessing.cpu_count()
        except Exception:
            available = 1
    n = max(1, math.floor(available * ratio))
    return n


def run_multiwfn():
    os.system("Multiwfn.exe < inputs.txt > Multi.log")


def generate_pdbs(file, working_dir, output_dir, cps):
    errors = []
    num_cpu = _detect_num_cpus(0.7)
    with open("inputs.txt", "w") as f:
        f.write(f"""
{file}
1000
10
{num_cpu}
100
2
2

2
1

0
2
{cps}
8
9
-20
-4
6
4
0
-5
6
0
-10
q
""")
    run_multiwfn()
    print("EXECUTION FINISHED")
    
    # Expected files
    has_paths = Path(working_dir, "paths.pdb").exists()
    has_cps_pdb = Path(working_dir, "CPs.pdb").exists()
    has_mol_pdb = Path(working_dir, f"{Path(file).stem}.pdb").exists()
    has_cps_txt = Path(working_dir, "CPs.txt").exists()

    if has_paths and has_cps_pdb and has_mol_pdb:
        shutil.move(Path(working_dir, "paths.pdb"), Path(output_dir, f"{Path(file).stem}_paths.pdb"))
        shutil.move(Path(working_dir, "CPs.pdb"), Path(output_dir, f"{Path(file).stem}_CPs.pdb"))
        shutil.move(Path(working_dir, f"{Path(file).stem}.pdb"), Path(output_dir, f"{Path(file).stem}.pdb"))
        shutil.move(Path(working_dir, f"{Path(file).stem}.xyz"), Path(output_dir, f"{Path(file).stem}.xyz"))
        
        # Move and rename the generated CPs.txt
        if has_cps_txt:
            shutil.move(Path(working_dir, "CPs.txt"), Path(output_dir, f"{Path(file).stem}_CPs.txt"))

        print(f"---> PDB and TXT files generated for structure, CPs and paths for {Path(file).name}")
    else:
        errors.append(f"{Path(file).stem}")
    if len(errors) > 0:
        print(f"!!! Errors found for {', '.join(errors)}")


def cp_properties(file, output_dir, cps, cp_indices, cps_txt_path):
    joined = f'\n7\n'.join(cp_indices)
    with open("inputs.txt", "w") as f:
        f.write(f"""
{file}
2
-4
5
{cps_txt_path}
0
7
{joined}
-10
q
""")
    os.system(f"Multiwfn.exe < inputs.txt > output.txt")
    with open('output.txt', 'r') as file2:
        lines = file2.readlines()
    start_indices = [i for i, line in enumerate(lines) if line.strip().startswith('Note: Unless otherwise specified, all units are')]
    end_indices = [i for i, line in enumerate(lines) if line.strip().startswith('eta index:')]
    if len(start_indices) != len(end_indices):
        print("!!! Error: Unequal number of 'start' and 'end' markers in CP results output")
        return
    i = 0
    for start_idx, end_idx in zip(start_indices, end_indices):
        extracted_lines = lines[start_idx:end_idx + 1]
        with open(Path(output_dir, f"{Path(file).stem}_cp_{cp_indices[i]}.txt"), "w") as output_file:
            output_file.writelines(extracted_lines)
        i += 1

# ==========================
# Main Script
# ==========================

working_dir = os.getcwd()

print("--- Directory Configuration ---")
wf_dir_input = input("Enter the path to the wavefunction folder (default: 'wf'): ").strip()
wf_dir_name = wf_dir_input.strip('"').strip("'") if wf_dir_input else 'wf'

output_dir_input = input("Enter the name or path of the output folder (default: 'output'): ").strip()
output_dir_name = output_dir_input.strip('"').strip("'") if output_dir_input else 'output'

cp_props_folder = 'cp_properties'
output_dir = Path(working_dir, output_dir_name)
cp_props_dir = Path(output_dir, cp_props_folder)
wf_dir = Path(working_dir, wf_dir_name)

print(f"\nBase working directory: {working_dir}")
print(f"Wavefunction directory: {wf_dir}")
print(f"Output directory: {output_dir}")

if not output_dir.exists(): os.makedirs(output_dir)
if not cp_props_dir.exists(): os.makedirs(cp_props_dir)

print("\nCritical point identification at the interface between fragments")

# Ask for execution mode
execution_mode = ""
while execution_mode not in ["1", "2"]:
    execution_mode = input(
        "Would you like to use Interactive or Batch mode for Ligand-Environment?\n"
        "[1] Interactive Mode (manual selection for each file)\n"
        "[2] Batch Mode (automatic ligand search across all files)\n> "
    ).strip()

ligand_xyz_path = ""
global_cp_types = {1, 2, 3, 4}
if execution_mode == "2":
    while True:
        ligand_xyz_path = input("Enter the path to the reference ligand .xyz file: ").strip()
        # strip quotes if user drags and drops
        ligand_xyz_path = ligand_xyz_path.strip('"').strip("'")
        if os.path.exists(ligand_xyz_path):
            break
        print("  ! File not found. Please try again.")
    global_cp_types = ime.read_cp_types_to_extract()

# ------------------------------ Section 1: Generate PDBs
cps = "2\n3\n4\n5"

# Collect all supported input files (.fchk, .wfx, .wfn)
wf_files = []
for ext in ["*.fchk", "*.wfx", "*.wfn"]:
    wf_files.extend(glob.glob(os.path.join(wf_dir, ext)))

base_names_wf = {Path(f).stem for f in wf_files}

print(f"Files found: {len(wf_files)}")

for file in wf_files:
    base = Path(file).stem
    print(f'Processing {base}...')

    existing = ime.check_existing_outputs(base, output_dir)
    # If all PDB and CSV already exist, skip this molecule
    if all(existing.values()):
        print("  - Existing PDB/CSV files detected; skipping PDB and CSV generation.")
        continue

    # Run local pdbs generation only if PDBs are missing
    if not (existing["pdb_structure"] and existing["pdb_cps"] and existing["pdb_paths"]):
        generate_pdbs(file, working_dir, output_dir, cps)
    else:
        print("  - PDBs already exist; skipping PDB generation.")

# ------------------------------ Section 2: Convert PDBs to CSV

all_pdb_files = glob.glob(os.path.join(output_dir, '*.pdb'))
# Structure PDBs
pdb_structure = [r for r in all_pdb_files
                 if not (Path(r).stem.endswith('_CPs') or Path(r).stem.endswith('_paths'))
                 and Path(r).stem in base_names_wf]
# CPs and paths PDBs
pdb_cps = [r for r in all_pdb_files if Path(r).stem.endswith('_CPs') and Path(r).stem[:-4] in base_names_wf]
pdb_paths = [r for r in all_pdb_files if Path(r).stem.endswith('_paths') and Path(r).stem[:-6] in base_names_wf]

for path in pdb_structure:
    dest = Path(output_dir, f"{Path(path).stem}.csv")
    if dest.exists():
        print(f"  - CSV already exists for {Path(path).stem}; skipping.")
    else:
        ime.csv_coordinates(working_dir, output_dir, path)

for path in pdb_cps:
    dest = Path(output_dir, f"{Path(path).stem}.csv")
    if dest.exists():
        print(f"  - CSV already exists for {Path(path).stem}; skipping.")
    else:
        ime.csv_cps(working_dir, output_dir, path)

for path in pdb_paths:
    dest = Path(output_dir, f"{Path(path).stem}.csv")
    if dest.exists():
        print(f"  - CSV already exists for {Path(path).stem}; skipping.")
    else:
        ime.csv_paths(working_dir, output_dir, path)

csv_files = glob.glob(os.path.join(output_dir, '*.csv'))
map_complex, map_cps, map_paths = {}, {}, {}
for csv_path in csv_files:
    stem = Path(csv_path).stem
    if stem.endswith('_CPs'):
        base = stem[:-4]
        if base in base_names_wf:
            map_cps[base] = csv_path
    elif stem.endswith('_paths'):
        base = stem[:-6]
        if base in base_names_wf:
            map_paths[base] = csv_path
    else:
        if stem in base_names_wf:
            map_complex[stem] = csv_path

if not map_complex:
    print("No structure files found to process.")
    raise SystemExit(0)

cp_indices_list, cp_vmd_indices_list, path_indices_list, processed_bases = [], [], [], []
overwrite = None

for base_name in sorted(map_complex.keys()):
    print(f"\nProcessing complex: {base_name}")
    cps_csv_path   = map_cps.get(base_name)
    paths_csv_path = map_paths.get(base_name)

    if cps_csv_path is None:
        print(f"  ! No CP CSV file found for {base_name}; skipping CP search.")
        cp_indices_list.append([]); cp_vmd_indices_list.append([]); path_indices_list.append([])
        continue

    complex_csv_path = map_complex[base_name]
    xyz_path = Path(output_dir, f"{base_name}.xyz")

    if execution_mode == "1":
        fragment_indices, cp_type, ignored_cps, visible_paths = ime.define_fragments_interactively(str(xyz_path))
    else:
        # Batch Mode
        fragment_indices = ime.define_fragments_batch(str(xyz_path), ligand_xyz_path)
        cp_type = global_cp_types
        ignored_cps = set()
        visible_paths = set()
    
    if not fragment_indices:
        print(f"  ! No fragments identified for {base_name}; skipping CP search.")
        cp_indices_list.append([]); cp_vmd_indices_list.append([]); path_indices_list.append([])
        continue

    coords, idx_nums = ime.load_csv_coordinates(complex_csv_path)
    # Load all CPs initially
    raw_cps_data = ime.load_csv_cps(cps_csv_path)
    
    # Filter out explicitly excluded CPs
    if ignored_cps:
        print(f"  - Excluding {len(ignored_cps)} manually marked CPs.")
        cps_data = [cp for cp in raw_cps_data if cp[0] not in ignored_cps]
    else:
        cps_data = raw_cps_data

    pair_to_cps, pair_to_vmd = ime._detect_interface_cps_no_cutoff(
        coords, idx_nums, cps_data, fragment_indices, valid_types=cp_type
    )

    interest_cps = set()
    interest_cps_vmd = set()
    
    # Print summary
    print(f"\nCP summary per interface for {base_name}:")
    for (fragA, fragB), cps_set in pair_to_cps.items():
        n_cps = len(cps_set)
        if n_cps > 0:
            print(f"  - {fragA} ↔ {fragB}: {n_cps} CPs")
            interest_cps.update([int(c) for c in cps_set])
            
            # VMD indices
            vmd_set = pair_to_vmd[(fragA, fragB)]
            interest_cps_vmd.update([int(c) for c in vmd_set])
        else:
             print(f"  - {fragA} ↔ {fragB}: 0 CPs  (—)")

    # Paths: 
    # In interactive mode, use those the user left visible.
    # In batch mode, extract automatically from interface CPs
    if execution_mode == "1":
        interest_path_vmd_indices = sorted([str(p) for p in visible_paths])
    else:
        interest_path_vmd_indices = []
        if paths_csv_path and interest_cps and not {2, 3, 4}.isdisjoint(cp_type):
            print("  - Searching for paths associated with the interface...")
            cp_types_map = {str(c[0]): int(c[1]) for c in cps_data}
            cps_for_path_search = [
                int(cp) for cp in interest_cps 
                if cp_types_map.get(str(cp)) != 1
            ]
            if cps_for_path_search:
                try:
                    paths = ime.path_indices_in_cps(
                        cps_csv_path, paths_csv_path, cps_for_path_search
                    )
                    interest_path_vmd_indices = sorted([str(p) for p in paths], key=lambda x: int(x))
                except Exception as e:
                    print(f"  ! Error getting paths for {base_name}: {e}")

    # Save to global lists
    cp_indices_list.append(list(interest_cps))
    cp_vmd_indices_list.append(list(interest_cps_vmd))
    path_indices_list.append(interest_path_vmd_indices)
    processed_bases.append(str(base_name))
    
    # Generate VMD Script
    if interest_cps_vmd:
        ime.load_structure_cps_paths(
            [f"{base_name}.pdb"],
            os.path.join(working_dir, 'load_structure_cps_paths.tcl'),
            [str(x) for x in interest_cps_vmd],
            interest_path_vmd_indices,
            str(output_dir)
        )

if cp_indices_list:
    for base_name, cps_list, cps_vmd_list, paths_list in zip(processed_bases, cp_indices_list, cp_vmd_indices_list, path_indices_list):
        cps_indices_str = " ".join([str(x) for x in sorted(set(cps_list), key=lambda x: int(x))] if cps_list else [])
        cps_vmd_str = " ".join([str(x) for x in sorted(set(cps_vmd_list), key=lambda x: int(x))] if cps_vmd_list else [])
        paths_str = " ".join([str(x) for x in sorted(set(paths_list), key=lambda x: int(x))] if paths_list else [])
        print(f"""
Complex {base_name}
        CPs found
        CP index list: {cps_indices_str if cps_indices_str else 'None'}
        CP index list for VMD (index - 1): {cps_vmd_str if cps_vmd_str else 'None'}
        Paths found
        Path index list for VMD (index): {paths_str if paths_str else 'None'}""")

overwrite = None
print("\nCalculating properties of extracted CPs")
if not cp_indices_list:
    print("No CPs identified; skipping property calculation.")
else:
    for base_name, cps_list in zip(processed_bases, cp_indices_list):
        if not cps_list:
            continue
        # Locate the wavefunction file
        wf_file_path = None

        for f in wf_files:
            if Path(f).stem == base_name:
                wf_file_path = Path(f)
                break
        
        if wf_file_path is None:
            print(f"No wavefunction file found for {base_name}; skipping its properties.")
            continue

        existing_pattern = str(Path(cp_props_dir, f"{base_name}_cp_*.txt"))
        try:
             cps_list_str = [str(x) for x in cps_list]
             specific_cps_txt_path = Path(output_dir, f"{base_name}_CPs.txt")
             cp_properties(str(wf_file_path), str(cp_props_dir), cps, cps_list_str, str(specific_cps_txt_path))
        except Exception as e:
            print(f"  ! Error calculating CP properties for {base_name}: {e}")

    complex_list = [f.stem for f in Path(output_dir).glob('*.tcl') if f.stem in base_names_wf]
    complex_list.sort()

    for complex_name in complex_list:
        cp_prop_txt_files = list(Path(cp_props_dir).glob('*_cp_*.txt'))
        matching_txt_files = [i for i in cp_prop_txt_files if i.stem.startswith(complex_name + '_cp_')]
        if not matching_txt_files:
            continue
        csv_prop_path = Path(cp_props_dir, f"{complex_name}_cp_properties.csv")
        
        try:
            ime.extract_cp_property_data(matching_txt_files, str(cp_props_dir))
        except Exception as e:
            print(f"  ! Error extracting property data for {complex_name}: {e}")
