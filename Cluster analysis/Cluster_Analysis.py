import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

def list_residues_and_clusters(u, cutoff_radius=10.0, cluster_size_threshold=5):
    # Initialize lists to store results for the selected frames
    all_results = []
    
    # Iterate over all frames in the trajectory
    for ts in u.trajectory:

        # Dictionary to store residue information
        residues_info = {}

        # Group atoms by residue
        for atom in u.atoms:
            residue_key = (atom.resid, atom.resname)
            if residue_key not in residues_info:
                residues_info[residue_key] = []
            residues_info[residue_key].append(atom)

        # Collect residues information
        res_atoms = []
        res_positions = []

        for (resid, resname), atom_list in residues_info.items():
            if resname in ["RES","SOL", "K", "Cl"]:
                # Get the position of the first atom in the residue for distance calculations
                res_atoms.append((resid, resname))
                res_positions.append(atom_list[0].position)

        # Calculate distances between all selected residues
        res_positions = np.array(res_positions)
        distance_matrix = distance_array(res_positions, res_positions)

        # Create adjacency matrix for residues within the cutoff radius
        adjacency_matrix = (distance_matrix < cutoff_radius).astype(int)

        # Find connected components (clusters)
        graph = csr_matrix(adjacency_matrix)
        n_components, labels = connected_components(csgraph=graph, directed=False)

        # Group residues by clusters
        clusters = {}
        for label, (resid, resname) in zip(labels, res_atoms):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append((resid, resname))

        # Filter clusters by size
        valid_clusters = [cluster for cluster in clusters.values() if len(cluster) >= cluster_size_threshold]

        # Track clustered residues
        clustered_residues = set()
        for cluster in valid_clusters:
            clustered_residues.update([resid for resid, resname in cluster])

        # Identify and count residues that do not participate in clustering
        non_clustered_residues = [resid for resid, resname in res_atoms if resid not in clustered_residues]
        num_non_clustered = len(non_clustered_residues)

        # Store results for the current frame
        frame_results = (ts.frame, len(valid_clusters), num_non_clustered)
        all_results.append(frame_results)

    return all_results

if __name__ == "__main__":
    # Replace with your actual file names
    tpr_file = '/home/12000_375mdrun10.tpr'
    xtc_file = '/home/12000_375mdrun10.xtc'

    # Load the topology and trajectory
    u = mda.Universe(tpr_file, xtc_file)

    # List residues information from the selected frames and find clusters
    all_results = list_residues_and_clusters(u)

    # Save results to a file
    output_file = '/home/ca1s_2000_375.xvg'
    with open(output_file, 'w') as f:
        f.write("frame_number\tnumber_of_cluster\tnumber_of_residues_notincluster\n")
        for frame_number, number_of_cluster, number_of_residues_notincluster in all_results:
            f.write(f"{frame_number}\t{number_of_cluster}\t{number_of_residues_notincluster}\n")

