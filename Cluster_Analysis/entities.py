import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix


def analyze_residues(u, cluster_size_threshold=5, small_fragment_radius=3.4):
    with open('result.xvg', 'w') as f:
        f.write("Frame Number\tNumber of Clusters\tNumber of Small Fragment Groups\t"
                "Number of Freely Roaming Fragments\tNumber of Freely Roaming RES Fragments\tNumber of Freely Roaming K Fragments\n")

        for ts in u.trajectory:
            if ts.frame == 0 or ts.frame % 1 == 0:
                box = np.append(ts.dimensions[:3], [90.0, 90.0, 90.0])

                residues_info = {}
                res_atoms = []
                res_positions = []

                for atom in u.atoms:
                    residue_key = (atom.resid, atom.resname)
                    if residue_key not in residues_info:
                        residues_info[residue_key] = []
                    residues_info[residue_key].append(atom)

                for (resid, resname), atom_list in residues_info.items():
                    if resname in ["RES", "SOL", "K", "Cl"]:
                        res_atoms.append((resid, resname))
                        res_positions.append(atom_list[0].position)

                if len(res_positions) == 0:
                    f.write(f"{ts.frame}\t0\t0\t0\t0\t0\n")
                    continue

                res_positions = np.array(res_positions)
                if res_positions.shape[0] > 0:
                    distance_matrix = distance_array(res_positions, res_positions, box=box)
                else:
                    f.write(f"{ts.frame}\t0\t0\t0\t0\t0\n")
                    continue

                has_near_neighbor = np.any(distance_matrix < small_fragment_radius, axis=1)
                adjacency_matrix = ((distance_matrix < 5.0) & has_near_neighbor[:, None]).astype(int)

                graph = csr_matrix(adjacency_matrix)
                n_components, labels = connected_components(csgraph=graph, directed=False)

                clusters = {}
                for label, (resid, resname) in zip(labels, res_atoms):
                    if label not in clusters:
                        clusters[label] = []
                    clusters[label].append((resid, resname))

                valid_clusters = [cluster for cluster in clusters.values() if len(cluster) >= cluster_size_threshold]
                number_of_clusters = len(valid_clusters)

                clustered_residues = set()
                for cluster in valid_clusters:
                    clustered_residues.update([resid for resid, _ in cluster])

                non_clustered_residues = [(resid, resname) for resid, resname in res_atoms if
                                          resid not in clustered_residues]
                if len(non_clustered_residues) == 0:
                    f.write(f"{ts.frame}\t{number_of_clusters}\t0\t0\t0\t0\n")
                    continue

                non_clustered_positions = [
                    res_positions[res_atoms.index((resid, resname))] for resid, resname in non_clustered_residues
                ]

                small_fragments = []
                freely_roaming_fragments = []
                freely_roaming_sol_fragments = 0
                freely_roaming_K_fragments = 0

                if len(non_clustered_positions) > 0:
                    for i, (resid, resname) in enumerate(non_clustered_residues):
                        neighbors = distance_array(
                            np.array(non_clustered_positions[i:i + 1]),
                            np.array(non_clustered_positions),
                            box=box
                        )[0]
                        neighbor_count = np.sum(neighbors < small_fragment_radius) - 1
                        if 1 <= neighbor_count <= 3:
                            small_fragments.append((resid, resname))
                        elif neighbor_count == 0:
                            freely_roaming_fragments.append((resid, resname))
                            if resname == "RES":
                                freely_roaming_sol_fragments += 1
                            elif resname == "K":
                                freely_roaming_K_fragments += 1

                if small_fragments:
                    small_fragment_positions = np.array([
                        non_clustered_positions[non_clustered_residues.index((resid, resname))]
                        for resid, resname in small_fragments
                    ])
                    if small_fragment_positions.shape[0] > 0:
                        small_fragment_distance_matrix = distance_array(small_fragment_positions,
                                                                        small_fragment_positions, box=box)
                        small_fragment_adjacency_matrix = (
                                    small_fragment_distance_matrix < small_fragment_radius).astype(int)
                        small_fragment_graph = csr_matrix(small_fragment_adjacency_matrix)
                        sf_components, sf_labels = connected_components(csgraph=small_fragment_graph, directed=False)
                        number_of_small_fragment_groups = sf_components
                    else:
                        number_of_small_fragment_groups = 0
                else:
                    number_of_small_fragment_groups = 0

                number_of_freely_roaming_fragments = len(freely_roaming_fragments)

                f.write(f"{ts.frame}\t{number_of_clusters}\t{number_of_small_fragment_groups}\t"
                        f"{number_of_freely_roaming_fragments}\t{freely_roaming_sol_fragments}\t{freely_roaming_K_fragments}\n")


if __name__ == "__main__":
    tpr_file = 'long375_10.tpr'
    xtc_file = 'long375_10.xtc'
    u = mda.Universe(tpr_file, xtc_file)
    analyze_residues(u)
