import MDAnalysis as mda
import numpy as np

def calculate_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def calculate_frame_counts(u, threshold=3.4, top_distances=6):
    # Dictionary to store oxygen and potassium atoms grouped by residue
    atoms_by_residue = {}

    # List to store the counts for the second frame
    frame_counts = []

    # Iterate over all frames
    for ts in u.trajectory:


        # Clear the atoms_by_residue dictionary in each iteration
        atoms_by_residue.clear()

        # Reset the count to zero for the frame
        count_condition_true = 0

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
        total_consecutive_res_count = 0
        consecutive_res_count = 0
        prev_res_number = None

        for residue_key_k in atoms_by_residue.keys():
            # Skip residues with no potassium atoms
            if 'K' not in [atom[0] for atom in atoms_by_residue[residue_key_k]]:
                continue

            k_atom_coord = np.array([atom[1] for atom in atoms_by_residue[residue_key_k] if atom[0] == 'K'][0])

            for residue_key_o in atoms_by_residue.keys():
                # Skip the same residue
                if residue_key_k == residue_key_o:
                    continue

                # Calculate distances and store in a list
                distances = []
                for atom_name, atom_coord in atoms_by_residue[residue_key_o]:
                    o_atom_coord = np.array(atom_coord)
                    distance = calculate_distance(k_atom_coord, o_atom_coord)

                    if distance <= threshold:
                        distances.append((distance, residue_key_o))

                if distances:
                    # Sort distances and get the top 6
                    sorted_distances = sorted(distances, key=lambda x: x[0])[:top_distances]

                    # Check for consecutive occurrences of the same RES number within the six smallest distances
                    res_number_counts = {}
                    for dist, res_key in sorted_distances:
                        current_res_number = res_key.split('RES')[0].strip()
                        if current_res_number in res_number_counts:
                            res_number_counts[current_res_number] += 1
                        else:
                            res_number_counts[current_res_number] = 1

                    # Print the count if it occurs six times consecutively
                    for res_number, count in res_number_counts.items():
                        if count >= 6:
                            consec_print_count += 1

                            # Track the total consecutive count
                            if res_number == prev_res_number:
                                consecutive_res_count += 1
                            else:
                                consecutive_res_count = 1
                            prev_res_number = res_number

                            if consecutive_res_count == 6:
                                total_consecutive_res_count += 1

        # Append the count for the second frame to the list
        frame_counts.append((ts.frame, consec_print_count))

    return frame_counts

if __name__ == "__main__":
    # Replace 'your_trajectory.tpr' and 'your_trajectory.xtc' with your actual file names
    tpr_file = '12000_375mdrun10.tpr'
    xtc_file = '12000_375mdrun10.xtc'

    # Load the topology and trajectory
    u = mda.Universe(tpr_file, xtc_file)

    # Process the trajectory and get the counts for the second frame
    frame_counts = calculate_frame_counts(u)

    # Save the counts for the second frame to a .xvg file
    output_file = 'cclast1s375_2000.xvg'
    np.savetxt(output_file, frame_counts, fmt='%d %d', header='Frame Count', comments='#')

