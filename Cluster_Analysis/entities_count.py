import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

def analyze_residues(u, cutoff_radius=10.0, cluster_size_threshold=5, small_fragment_radius=5.0):
    output_file = 'entites.xvg'

    with open(output_file, 'w') as f:
        header = ("Frame Number\tFreely Roaming Fragments\tFreely Roaming K\tSmall Fragment Groups\t"
                  "Number of Agglomerates\n")
        f.write(header)
        print(header.strip())

        for ts in u.trajectory:
            if ts.frame % 1 != 0:
                continue
            box = np.append(ts.dimensions[:3], [90.0, 90.0, 90.0])

            res_atoms = []
            res_positions = []

            for residue in u.residues:
                if residue.resname in ["RES", "SOL"]:
                    com = residue.atoms.center_of_mass()
                    res_atoms.append((residue.resid, residue.resname))
                    res_positions.append(com)
                elif residue.resname in ["K", "Cl"]:
                    res_atoms.append((residue.resid, residue.resname))
                    res_positions.append(residue.atoms[0].position)

            if len(res_positions) == 0:
                result_line = f"{ts.frame}\t0\t0\t0\t0\n"
                f.write(result_line)
                print(result_line.strip())
                continue

            res_positions = np.array(res_positions)


            freely_roaming_fragments = []
            freely_roaming_K = 0
            filtered_atoms = []
            filtered_positions = []

            distance_matrix = distance_array(res_positions, res_positions, box=box)

            for i, (resid, resname) in enumerate(res_atoms):
                neighbors = distance_matrix[i] < small_fragment_radius
                if np.sum(neighbors) == 1:  # Only self
                    freely_roaming_fragments.append((resid, resname))
                    if resname == "K":
                        freely_roaming_K += 1
                else:
                    filtered_atoms.append((resid, resname))
                    filtered_positions.append(res_positions[i])

            res_atoms, res_positions = filtered_atoms, np.array(filtered_positions)

            if len(res_positions) == 0:
                result_line = f"{ts.frame}\t{len(freely_roaming_fragments)}\t{freely_roaming_K}\t0\t0\n"
                f.write(result_line)
                print(result_line.strip())
                continue


            small_fragments = []
            adjacency_matrix = (distance_array(res_positions, res_positions, box=box) < small_fragment_radius).astype(int)

            if adjacency_matrix.size > 0:
                graph = csr_matrix(adjacency_matrix)
                n_components, labels = connected_components(csgraph=graph, directed=False)
                small_fragment_groups = {i: [] for i in range(n_components)}

                for label, (resid, resname) in zip(labels, res_atoms):
                    small_fragment_groups[label].append((resid, resname))

                small_fragments = [group for group in small_fragment_groups.values() if 1 <= len(group) <= 3]

            small_fragment_residues = set(resid for group in small_fragments for resid, _ in group)
            filtered_atoms = []
            filtered_positions = []

            for i, (resid, resname) in enumerate(res_atoms):
                if resid not in small_fragment_residues:
                    filtered_atoms.append((resid, resname))
                    filtered_positions.append(res_positions[i])

            res_atoms, res_positions = filtered_atoms, np.array(filtered_positions)

            if len(res_positions) == 0:
                result_line = f"{ts.frame}\t{len(freely_roaming_fragments)}\t{freely_roaming_K}\t{len(small_fragments)}\t0\n"
                f.write(result_line)
                print(result_line.strip())
                continue


            number_of_clusters = 0
            if len(res_positions) > 0:
                adjacency_matrix = (distance_array(res_positions, res_positions, box=box) < cutoff_radius).astype(int)
                graph = csr_matrix(adjacency_matrix)
                n_components, labels = connected_components(csgraph=graph, directed=False)
                clusters = {i: [] for i in range(n_components)}

                for label, (resid, resname) in zip(labels, res_atoms):
                    clusters[label].append((resid, resname))

                valid_clusters = [cluster for cluster in clusters.values() if len(cluster) >= cluster_size_threshold]
                number_of_clusters = len(valid_clusters)

            result_line = f"{ts.frame}\t{len(freely_roaming_fragments)}\t{freely_roaming_K}\t{len(small_fragments)}\t{number_of_clusters}\n"
            f.write(result_line)
            print(result_line.strip())

if __name__ == "__main__":
    tpr_file = 'xxxx.tpr'
    xtc_file = 'xxxx.xtc'
    u = mda.Universe(tpr_file, xtc_file)
    analyze_residues(u)
