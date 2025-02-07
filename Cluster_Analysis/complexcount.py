import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array

def calculate_distance(coord1, coord2, box):
    """Calculate distance using periodic boundary conditions (PBC)."""
    dist = distance_array(coord1[np.newaxis, :], coord2[np.newaxis, :], box=box)
    return dist[0, 0]

def calculate_frame_counts(u, threshold=3.4, top_distances=6, frame_range=(0, 10000)):
    # Dictionary to store oxygen and potassium atoms grouped by residue
    atoms_by_residue = {}

    # List to store detailed results for each frame
    frame_counts = []

    # Iterate over frames within the specified range
    for ts in u.trajectory:
        if not (frame_range[0] <= ts.frame <= frame_range[1]):
            continue

        # Only process frames divisible by 5
        if ts.frame % 5 != 0:
            continue

        # Clear the atoms_by_residue dictionary in each iteration
        atoms_by_residue.clear()

        # Reset the count to zero for each frame
        count_condition_true = 0

        # Get the simulation box dimensions and apply PBC
        box = np.append(ts.dimensions[:3], [90.0, 90.0, 90.0])

        # Select oxygen atoms with names starting with 'O'
        oxygen_atoms = u.select_atoms('name O* and not name OW*')

        # Group oxygen atoms by residue
        for atom in oxygen_atoms:
            residue_key = f'{atom.resid}RES'
            if 43929 <= atom.resid <= 44428 and residue_key not in atoms_by_residue:
                atoms_by_residue[residue_key] = []
            atoms_by_residue[residue_key].append((atom.name, atom.position))

        # Select potassium atoms with name 'K'
        potassium_atoms = u.select_atoms('name K')

        # Group potassium atoms by residue
        for atom in potassium_atoms:
            residue_key = f'{atom.resid}K'  # Assuming potassium residue numbers are not modified
            if 42929 <= atom.resid <= 43927 and residue_key not in atoms_by_residue:
                atoms_by_residue[residue_key] = []
            atoms_by_residue[residue_key].append((atom.name, atom.position))

        # Calculate distances and average distances
        consec_print_count = 0
        qualifying_k_residues = []  # To store potassium residues satisfying the condition

        for residue_key_k in atoms_by_residue.keys():
            if 'K' not in [atom[0] for atom in atoms_by_residue[residue_key_k]]:
                continue

            k_atom_coord = np.array([atom[1] for atom in atoms_by_residue[residue_key_k] if atom[0] == 'K'][0])

            for residue_key_o in atoms_by_residue.keys():
                if residue_key_k == residue_key_o:
                    continue

                distances = []
                for atom_name, atom_coord in atoms_by_residue[residue_key_o]:
                    o_atom_coord = np.array(atom_coord)
                    distance = calculate_distance(k_atom_coord, o_atom_coord, box)

                    if distance <= threshold:
                        distances.append((distance, residue_key_o))

                if distances:
                    sorted_distances = sorted(distances, key=lambda x: x[0])[:top_distances]

                    res_number_counts = {}
                    for dist, res_key in sorted_distances:
                        current_res_number = res_key.split('RES')[0].strip()
                        res_number_counts[current_res_number] = res_number_counts.get(current_res_number, 0) + 1

                    for res_number, count in res_number_counts.items():
                        if count >= 6:
                            consec_print_count += 1
                            qualifying_k_residues.append(residue_key_k)

        hydrated_count = 0
        non_hydrated_count = 0
        hydrated_residues = []
        non_hydrated_residues = []

        sol_atoms = u.select_atoms('resname SOL')

        for residue_key_k in qualifying_k_residues:
            k_atom_coord = np.array([atom[1] for atom in atoms_by_residue[residue_key_k] if atom[0] == 'K'][0])

            hydrated = False
            for sol_atom in sol_atoms:
                sol_coord = sol_atom.position
                distance = calculate_distance(k_atom_coord, sol_coord, box)
                if distance <= 3.4:
                    hydrated = True
                    break

            if hydrated:
                hydrated_count += 1
                hydrated_residues.append(residue_key_k)
            else:
                non_hydrated_count += 1
                non_hydrated_residues.append(residue_key_k)

        frame_counts.append({
            "frame": ts.frame,
            "consec_count": consec_print_count,
            "hydrated_count": hydrated_count,
            "non_hydrated_count": non_hydrated_count,
            "hydrated_residues": hydrated_residues,
            "non_hydrated_residues": non_hydrated_residues,
        })

    return frame_counts

if __name__ == "__main__":
    tpr_file = '12.tpr'
    xtc_file = '12.xtc'

    # Load the topology and trajectory
    u = mda.Universe(tpr_file, xtc_file)

    # Process the trajectory and get the counts for each frame in the range 0-10,000
    frame_counts = calculate_frame_counts(u, frame_range=(0, 10000))

    # Save the counts for each frame to a .xvg file
    output_file = 'result.xvg'
    with open(output_file, 'w') as f:
        f.write("# Frame Consecutive_Count Hydrated_Count Non_Hydrated_Count\n")
        for entry in frame_counts:
            frame = entry["frame"]
            consec_count = entry["consec_count"]
            hydrated_count = entry["hydrated_count"]
            non_hydrated_count = entry["non_hydrated_count"]
            f.write(f"{frame} {consec_count} {hydrated_count} {non_hydrated_count}\n")
