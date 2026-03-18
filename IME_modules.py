"""
Module functions for Multiwfn and PDB/CSV/CP/path handling.
"""
import re
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
import os, math
import multiprocessing
import glob
import itertools
from collections import defaultdict
import FragmentFinder as Ff


def csv_coordinates(working_dir, output_dir, filepath):
    try:
        from ase.io import read
        atoms = read(filepath)
        
        data = []
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        for i, (pos, sym) in enumerate(zip(positions, symbols)):
            data.append([i + 1, sym, pos[0], pos[1], pos[2]])

        df = pd.DataFrame(data, columns=['number', 'atom', 'x', 'y', 'z'])

        out = Path(output_dir) / f"{Path(filepath).stem}.csv"
        df.to_csv(out, index=False)
        if out.exists():
            print(f"---> CSV file generated for coordinates of {Path(filepath).stem} (via ASE)")

    except Exception as e:
        print(f"!!! Critical error reading PDB with ASE: {e}. Make sure 'ase' is installed correctly.")


def csv_cps(working_dir, output_dir, filepath):
    data = []
    with open(Path(filepath), 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                CP = line[6:11].strip()
                elem = line[76:78].strip() if len(line) >= 78 else ''
                # Map element to type
                cp_type = '1' if elem == 'C' else '2' if elem == 'N' else '3' if elem == 'O' else '4' if elem == 'F' else 'none'
                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()
                data.append([CP, cp_type, x, y, z])

    df = pd.DataFrame(data, columns=['CP', 'tipo', 'x', 'y', 'z'])
    df['CP'] = pd.to_numeric(df['CP'], errors='coerce')
    for c in ['x','y','z']:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    out = Path(output_dir) / f"{Path(filepath).stem}.csv"
    df.to_csv(out, index=False)
    if out.exists():
        print(f"---> CSV file generated for coordinates of {Path(filepath).name}")
    else:
        print(f'!!! CSV write error for {Path(filepath).stem}')


def csv_paths(working_dir, output_dir, filepath):
    data = []
    with open(Path(filepath), 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                number = line[22:26].strip()
                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()
                data.append([number, x, y, z])

    df = pd.DataFrame(data, columns=['number', 'x', 'y', 'z'])
    df['number'] = df['number'].astype(str)
    for c in ['x','y','z']:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    out = Path(output_dir) / f"{Path(filepath).stem}.csv"
    df.to_csv(out, index=False)
    if out.exists():
        print(f"---> CSV file generated for coordinates of {Path(filepath).name}")
    else:
        print(f'!!! CSV write error for {Path(filepath).stem}')


def path_indices_in_cps(CPs_csv, paths_csv, cp_list):
    # Cutoff value to find points close to CPs
    cutoff = 0.02

    CPs_df = pd.read_csv(CPs_csv)
    paths_df = pd.read_csv(paths_csv)

    path_list = []

    # For each CP, calculate the distance to each path point
    for cp_row_idx, cp_row in CPs_df.iterrows():
        if int(cp_row['CP']) in cp_list:
            for path_row_idx, path_row in paths_df.iterrows():
                distance = np.sqrt((cp_row['x'] - path_row['x'])**2 +
                                    (cp_row['y'] - path_row['y'])**2 +
                                    (cp_row['z'] - path_row['z'])**2)

                if distance <= cutoff:
                    path_list.append(str(int(path_row['number'])))

    return path_list


def load_structure_cps_paths(file_list, tcl_template_path, interest_cps, interest_paths, output_dir):
    if not os.path.isfile(tcl_template_path):
        print(f"!!! TCL template not found at '{tcl_template_path}'. Scripts will not be generated.")
        return
    with open(tcl_template_path, "r") as f_tpl:
        tcl_content = f_tpl.read()
    cps_vmd = " ".join(interest_cps) if interest_cps else "none"
    
    if interest_paths:
         paths_selection_str = f"resid {' '.join(interest_paths)}"
         show_paths_rep = "1" # Show
    else:
         paths_selection_str = "none"
         show_paths_rep = "0" # Hide

    for filename in file_list:
        pdb_name = Path(filename).name
        stem = Path(filename).stem
        cp_name = f"{stem}_CPs.pdb"
        paths_name = f"{stem}_paths.pdb"
        mod = tcl_content
        mod = mod.replace("[structure file]", pdb_name)\
                 .replace("[CPs file]", cp_name)\
                 .replace("[paths file]", paths_name)\
                 .replace("[interest_cps]", cps_vmd)\
                 .replace("[selection_paths]", paths_selection_str)\
                 .replace("[show_paths_rep]", show_paths_rep)
        new_name = f"{stem}.tcl"
        new_tcl_path = os.path.join(output_dir, new_name)
        with open(new_tcl_path, "w") as out:
            out.write(mod)
        print(f"---> File {Path(new_tcl_path).name} created in {Path(output_dir).name} successfully.")


def extract_cp_property_data(files, output_dir):

    properties = [
        "CP type", "Density of all electrons", "Density of Alpha electrons",
        "Density of Beta electrons", "Spin density of electrons",
        "Lagrangian kinetic energy G(r)", "Hamiltonian kinetic energy K(r)",
        "Potential energy density V(r)", "Energy density E(r) or H(r)",
        "Laplacian of electron density", "Electron localization function (ELF)",
        "Localized orbital locator (LOL)", "Local information entropy",
        "Interaction region indicator (IRI)", "Reduced density gradient (RDG)",
        "Reduced density gradient with promolecular approximation",
        "Sign(lambda2)*rho", "Sign(lambda2)*rho with promolecular approximation",
        "Wavefunction value for orbital       1", "Average local ionization energy (ALIE)",
        "van der Waals potential (probe atom: C )", "Delta-g (under promolecular approximation)",
        "Delta-g (under Hirshfeld partition)", "ESP from nuclear charges",
        "ESP from electrons", "Total ESP"
    ]

    num_re = re.compile(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?')

    def _first_numeric(s):
        m = num_re.search(s)
        if not m:
            return np.nan
        token = m.group(0).replace('d', 'E').replace('D', 'E')
        try:
            return float(token)
        except Exception:
            return np.nan

    data = {prop: [] for prop in properties}
    cp_numbers = []

    for filepath in files:
        filepath = Path(filepath)
        cpnum = filepath.stem.split('_')[-1].split('.')[0]
        cp_numbers.append(cpnum)

        try:
            lines = filepath.read_text(encoding='utf-8', errors='ignore').splitlines()
        except Exception:
            with open(filepath, 'r', errors='ignore') as f:
                lines = f.readlines()

        # For each property, find its line and extract the value
        for prop in properties:
            value = None
            prefix = prop
            for line in lines:
                s = line.strip()
                if s.startswith(prefix):
                    rhs = s.split(':', 1)[1].strip() if ':' in s else s[len(prefix):].strip()
                    if prop == "CP type":
                        if '(' in rhs and ')' in rhs:
                            inside = rhs[rhs.find('(')+1 : rhs.find(')')]
                        else:
                            inside = rhs
                        value = inside.replace(' ', '')
                    else:
                        value = _first_numeric(rhs)
                    break
            if value is None:
                data[prop].append(np.nan if prop != "CP type" else "")
            else:
                data[prop].append(value)

    df = pd.DataFrame(data, index=[str(i) for i in cp_numbers])
    df.index.name = 'ID'

    # Sort rows by CP number (index)
    try:
        df['sort_val'] = df.index.astype(int)
        df = df.sort_values('sort_val')
        df = df.drop('sort_val', axis=1)
    except Exception:
        pass

    # Write CSV
    name = Path(files[0]).stem.split('_cp_')[0]
    out_path = Path(output_dir) / f'{name}_cp_properties.csv'
    df.to_csv(out_path, header=True)
    print(f"CP properties CSV file created successfully in {out_path.parent.name}")

    return df, out_path


def define_fragments_interactively(xyz_path, valid_cp_types=None):
    """
    Launches the unified graphical interface to define fragments and visualize CPs.
    """
    path_obj = Path(xyz_path)
    name = path_obj.stem
    folder = path_obj.parent
    
    # Load CP data
    cps_path = folder / f"{name}_CPs.csv"
    cp_data = []
    if cps_path.exists():
        try:
            cp_data = load_csv_cps(str(cps_path))
        except Exception as e:
            print(f"Warning: Could not load CPs for visualization: {e}")

    # Load path data
    paths_path = folder / f"{name}_paths.csv"
    path_data = {}
    if paths_path.exists():
        try:
            df_p = pd.read_csv(paths_path)
            groups = df_p.groupby('number')
            for pid, group in groups:
                coords = group[['x','y','z']].to_numpy()
                path_data[int(pid)] = coords
        except Exception as e:
             print(f"Warning: Could not load paths for visualization: {e}")

    # Analysis callback
    # This callback runs inside FragmentFinder when the user presses 'f'
    def analysis_callback(coords, idx_nums, raw_cps, frag_indices_dict, active_types):
        return _detect_interface_cps_no_cutoff(
            coords, idx_nums, raw_cps, frag_indices_dict,
            valid_types=active_types
        )

    # Interactive session
    print(f"\nStarting interactive selector for {name}...")
    try:
        fragment_indices, selected_types, ignored_cps, visible_paths = Ff.run_interactive_session(xyz_path, cp_data, path_data, analysis_callback)
    except AttributeError:
        print("Error: FragmentFinder not updated. Using fallback method.")
        return {}, {1,2,3,4}, set(), set()
    except ValueError:
        print("Warning: Old return signature detected.")
        return {}, {1,2,3,4}, set(), set()

    return fragment_indices, selected_types, ignored_cps, visible_paths


def define_fragments_batch(complex_xyz: str, ligand_xyz: str):
    """
    Automatically searches for the reference ligand within the complex using graph isomorphism.
    Assigns the ligand to Fragment 1 and the remaining atoms to Fragment 2.
    """
    try:
        from ase.io import read
        complex_mol = read(complex_xyz)
        ligand_mol = read(ligand_xyz)
    except Exception as e:
        print(f"  ! Error reading XYZ files: {e}")
        return {}
        
    complex_matrix = Ff.calculate_connectivity_matrix(complex_mol)
    complex_symbols = complex_mol.get_chemical_symbols()
    
    ligand_matrix = Ff.calculate_connectivity_matrix(ligand_mol)
    ligand_symbols = ligand_mol.get_chemical_symbols()
    
    print(f"  - Searching for ligand ({len(ligand_symbols)} atoms) in complex ({len(complex_symbols)} atoms)...")
    matches = Ff.match_fragment(complex_matrix, ligand_matrix, complex_symbols, ligand_symbols)
    
    if not matches:
        print("  ! No matching topology found for the ligand in this complex.")
        return {}
    match_indices_0 = matches[0]
    
    match_set_0 = set(match_indices_0)
    all_set_0 = set(range(len(complex_symbols)))
    remaining_set_0 = all_set_0 - match_set_0
    
    frag_1 = {i + 1 for i in match_set_0}
    frag_2 = {i + 1 for i in remaining_set_0}
    
    fragment_indices = {1: frag_1, 2: frag_2}
    print(f"  - Ligand successfully mapped to Fragment 1.")
    print(f"  - Environment assigned to Fragment 2 ({len(frag_2)} atoms).")
    
    return fragment_indices


def _detect_interface_cps_no_cutoff(coords, idx_atom_numbers, cps_data,
                                     fragment_indices, valid_types,
                                     k_neighbors=4, eps=1e-12):
    """
    Assigns CPs to interfaces (A,B) without using a distance cutoff.

    Based on:
      1) Local neighborhood: takes the k closest atoms to P. If all belong to
         the same fragment => P is "inside" that fragment => discard.
      2) For each pair of fragments present in the neighborhood, finds the pair
         that MINIMIZES the distance from P to the segment a-b.
      3) Accepts if:
           - the projection falls within the segment
           - the vectors P->a and P->b are approximately opposite.
    """
    valid_types = set(valid_types) if valid_types else None

    # Fast lookup maps
    number_to_pos = {int(num): i for i, num in enumerate(idx_atom_numbers)}
    atom_to_frag  = {}
    for frag, atoms in fragment_indices.items():
        for n in atoms:
            atom_to_frag[int(n)] = frag

    # Positions and coords per fragment
    frag_pos, frag_coords = {}, {}
    for frag, atoms in fragment_indices.items():
        pos = [number_to_pos[n] for n in atoms if n in number_to_pos]
        if pos:
            frag_pos[frag]    = np.asarray(pos, dtype=int)
            frag_coords[frag] = coords[frag_pos[frag], :]

    if len(frag_coords) < 2:
        return defaultdict(set), defaultdict(set)

    all_atom_frags = np.array([atom_to_frag.get(int(n), None) for n in idx_atom_numbers], dtype=object)

    pair_to_cps = defaultdict(set)
    pair_to_vmd = defaultdict(set)

    for cp_id, cp_type, cp_xyz in cps_data:
        if valid_types and (cp_type not in valid_types):
            continue

        P = np.asarray(cp_xyz, dtype=float)
        if cp_type == 1:
            d2_all = np.sum((coords - P) ** 2, axis=1)
            nearest_idx = np.argmin(d2_all)
            nearest_dist = np.sqrt(d2_all[nearest_idx])
            
            # Threshold for "coincident" (0.1 Angstrom for nuclear CPs)
            if nearest_dist < 0.1:
                real_atom_num = idx_atom_numbers[nearest_idx]
                owner_frag = atom_to_frag.get(int(real_atom_num))
                
                # If we found an owner, this CP is strictly internal
                if owner_frag is not None:
                     continue

        if cp_type == 2:
            k_dynamic = 2
        elif cp_type == 3:
            k_dynamic = 3
        else:
            k_dynamic = 4 

        # Local neighborhood of P
        d2_all = np.sum((coords - P) ** 2, axis=1)
        k = min(k_dynamic, d2_all.shape[0])
        nn_idx = np.argpartition(d2_all, k-1)[:k]
        nn_frags = {f for i in nn_idx if (f := all_atom_frags[i]) is not None}

        # If only 1 fragment in neighborhood => internal => discard
        if len(nn_frags) < 2:
            continue

        # Evaluate only pairs present in the neighborhood
        for A, B in itertools.combinations(sorted(nn_frags), 2):
            Acoords = frag_coords.get(A); Bcoords = frag_coords.get(B)
            if Acoords is None or Bcoords is None:
                continue

            # Distance from P to the best segment a-b
            ab = Bcoords[None, :, :] - Acoords[:, None, :]
            ap = P[None, None, :]  - Acoords[:, None, :]

            denom = np.sum(ab * ab, axis=2) + eps
            t     = np.sum(ap * ab, axis=2) / denom 
            # Closest point on the segment
            t_clip = np.clip(t, 0.0, 1.0)
            proj   = Acoords[:, None, :] + t_clip[..., None] * ab
            dist2  = np.sum((proj - P) ** 2, axis=2)

            # Best pair by minimum distance to segment
            ia, ib = np.unravel_index(np.argmin(dist2), dist2.shape)
            t_best = float(t[ia, ib])

            # Geometric test: within segment + angular opposition
            if 0.0 <= t_best <= 1.0:
                a = Acoords[ia]; b = Bcoords[ib]
                v1 = a - P; v2 = b - P
                n1 = np.linalg.norm(v1); n2 = np.linalg.norm(v2)
                if n1 > eps and n2 > eps:
                    cosang = float(np.dot(v1, v2) / (n1 * n2))
                    if cosang <= 0.0: 
                        pair = (A, B) if A < B else (B, A)
                        pair_to_cps[pair].add(str(cp_id))
                        pair_to_vmd[pair].add(str(cp_id - 1))

    return pair_to_cps, pair_to_vmd


def load_csv_coordinates(csv_path):
    """
    Loads a CSV with columns: number, atom, x, y, z
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    df = pd.read_csv(path, sep=None, engine="python",
                     encoding="utf-8-sig", comment="#")

    if df.shape[1] == 1:
        df = pd.read_csv(path, delim_whitespace=True, engine="python",
                         encoding="utf-8-sig", comment="#")

    df.columns = [str(c).strip().lower() for c in df.columns]

    colmap = {}
    if "number" not in df.columns:
        for c in ("numero", "num", "index", "idx", "id"):
            if c in df.columns:
                colmap[c] = "number"; break
    if "atom" not in df.columns:
        for c in ("atomo", "symbol", "elemento", "elem"):
            if c in df.columns:
                colmap[c] = "atom"; break
    df = df.rename(columns=colmap)
    req = {"number", "x", "y", "z"}
    missing = req.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {path.name}: {sorted(missing)}")

    # Cleanup and conversion
    for c in ("x", "y", "z"):
        df[c] = (df[c].astype(str).str.replace(",", ".", regex=False))
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["number"] = pd.to_numeric(df["number"], errors="coerce")

    # Remove badly parsed rows
    df = df.dropna(subset=["number", "x", "y", "z"])
    df["number"] = df["number"].astype(int)

    # Sort by 'number'
    df = df.sort_values("number")

    coords = df[["x", "y", "z"]].to_numpy(dtype=float)
    idx_nums = df["number"].to_numpy(dtype=int)

    return coords, idx_nums


def load_csv_cps(cps_csv_path: str):
    df = pd.read_csv(cps_csv_path)
    cps = []
    for _, r in df.iterrows():
        try:
            cp = int(r['CP']); t = int(r['tipo'])
            xyz = np.array([float(r['x']), float(r['y']), float(r['z'])], dtype=float)
            cps.append((cp, t, xyz))
        except Exception:
            continue
    return cps


def read_cp_types_to_extract():
    """Prompts the user for the CP types to extract"""
    valid = {1, 2, 3, 4}
    while True:
        txt = input(
            "\nSelect the Critical Point (CP) types to extract:\n"
            "  [1] (3,-3) Nuclear  - Nucleus (VMD Color: Silver)\n"
            "  [2] (3,-1) Bond     - Bond    (VMD Color: Lime)\n"
            "  [3] (3,+1) Ring     - Ring    (VMD Color: Purple)\n"
            "  [4] (3,+3) Cage     - Cage    (VMD Color: Yellow)\n\n"
            "Enter the numbers separated by commas (e.g. 2,3) (Include 2 to view Paths): "
        ).strip()
        try:
            values = {int(t.strip()) for t in txt.split(",") if t.strip()}
        except Exception:
            print("  ! Invalid input. Enter numbers separated by commas in the range 1-4.")
            continue
        if not values or not values.issubset(valid):
            print("  ! Only values 1,2,3,4 are allowed.")
            continue
        return values

def check_existing_outputs(base, output_dir):
    """Checks if the main PDB/CSV files exist for a molecule `base`."""
    expect = {
        "pdb_structure": output_dir / f"{base}.pdb",
        "pdb_cps":       output_dir / f"{base}_CPs.pdb",
        "pdb_paths":     output_dir / f"{base}_paths.pdb",
        "csv_structure": output_dir / f"{base}.csv",
        "csv_cps":       output_dir / f"{base}_CPs.csv",
        "csv_paths":     output_dir / f"{base}_paths.csv",
    }
    existing = {k: p.exists() for k, p in expect.items()}
    return existing
