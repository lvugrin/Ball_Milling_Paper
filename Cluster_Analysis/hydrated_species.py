import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix


def analyze_residues(u):
    output_file = 'result.xvg'

    with open(output_file, 'w') as f:
        header = ("Frame Number\tHydrated Crowns\t"
                  "Crowns_with_Crowns_withwater\tCrowns_with_Crowns_nowater\tHydrated_K_alone\t"
                  "HydratedK_kcl\tK_with_Cl\n")
        f.write(header)
        print(header.strip())

        for ts in u.trajectory:
            if ts.frame % 5 != 0:
                continue
            box = np.append(ts.dimensions[:3], [90.0, 90.0, 90.0])

            residues_info = {}
            res_atoms = []
            res_positions = []

            # Group atoms by residue and calculate COM where applicable
            for residue in u.residues:
                residue_key = (residue.resid, residue.resname)
                if residue.resname in ["RES", "SOL"]:
                    com = residue.atoms.center_of_mass()
                    res_atoms.append((residue.resid, residue.resname))
                    res_positions.append(com)
                elif residue.resname in ["K", "Cl"]:
                    res_atoms.append((residue.resid, residue.resname))
                    res_positions.append(residue.atoms[0].position)

            if len(res_positions) == 0:
                result_line = f"{ts.frame}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
                f.write(result_line)
                print(result_line.strip())
                continue

            res_positions = np.array(res_positions)
            if res_positions.shape[0] > 0:
                distance_matrix = distance_array(res_positions, res_positions, box=box)
            else:
                result_line = f"{ts.frame}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
                f.write(result_line)
                print(result_line.strip())
                continue

            RES_with_SOL_within_4_no_K_no_RES = 0
            RES_with_RES_or_Cl_no_SOL_or_K = 0
            K_with_SOL_within_4_no_RES = 0
            K_with_SOL_within_4_no_RES_no_otherK = 0
            K_with_K_or_Cl_no_SOL_or_RES = 0
            RES_with_SOL_within_4_no_K_yes_RES = 0
            has_SOL_within_4_count_RES = 0
            has_SOL_within_4_count_K = 0
            RES_with_RES_or_Cl_yes_SOL_or_K = 0
            K_with_SOL_within_4_no_RES_has_otherK_Cl = 0


            for i, (resid, resname) in enumerate(res_atoms):
                if resname == "RES":
                    # Calculate distances to various residue types
                    distances_to_SOL = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "SOL"
                    ]
                    distances_to_K = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "K"
                    ]
                    distances_to_otherRES = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "RES"
                    ]
                    distances_to_Cl = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "Cl"
                    ]

                    # Count RES_with_SOL_within_4_no_K
                    has_SOL_within_4 = False
                    for d in distances_to_SOL:
                        if d <= 3.4:
                            has_SOL_within_4 = True
                            break

                    no_K_within_4 = True
                    for d in distances_to_K:
                        if d <= 3.0:
                            no_K_within_4 = False
                            break

                    no_RES_within_4 = True
                    for d in distances_to_otherRES:
                        if 0 < d <= 4.0:
                            no_RES_within_4 = False
                            break

                    if has_SOL_within_4 and no_K_within_4 and no_RES_within_4:
                        RES_with_SOL_within_4_no_K_no_RES += 1

                    # Count RES_with_RES_or_Cl_no_SOL_or_K
                    no_SOL_within_4 = True
                    for d in distances_to_SOL:
                        if d <= 4.0:
                            no_SOL_within_4 = False
                            break

                    no_K_within_4 = True
                    for d in distances_to_K:
                        if d <= 3.0:
                            no_K_within_4 = False
                            break

                    has_RES_within_4 = False
                    for d in distances_to_otherRES:
                        if 0 < d <= 5.0:
                            has_RES_within_4 = True
                            break

                    if no_SOL_within_4 and no_K_within_4 and has_RES_within_4:
                        RES_with_RES_or_Cl_no_SOL_or_K += 1

                    has_SOL_within_6 = False
                    for d in distances_to_SOL:
                        if d <= 5.0:
                            has_SOL_within_6 = True
                            break
                    no_K_within_4 = True
                    for d in distances_to_K:
                        if d <= 3.0:
                            no_K_within_4 = False
                            break

                    has_RES_within_4 = False
                    for d in distances_to_otherRES:
                        if 0 < d <= 5.0:
                            has_RES_within_4 = True
                            break

                    if has_SOL_within_6 and no_K_within_4 and has_RES_within_4:
                        RES_with_RES_or_Cl_yes_SOL_or_K += 1

            for i, (resid, resname) in enumerate(res_atoms):
                if resname == "K":
                    # Calculate distances to various residue types
                    distances_to_SOL = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "SOL"
                    ]
                    distances_to_RES = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "RES"
                    ]
                    distances_to_otherK = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "K"
                    ]
                    distances_to_Cl = [
                        distance_matrix[i, j]
                        for j, (_, rn) in enumerate(res_atoms) if rn == "Cl"
                    ]

                    # Count K_with_SOL_within_4_no_RES
                    has_SOL_within_3_4 = False
                    for d in distances_to_SOL:
                        if d <= 3.4:
                            has_SOL_within_3_4 = True
                            break

                    no_RES_within_3_4 = True
                    for d in distances_to_RES:
                        if d <= 3.0:
                            no_RES_within_3_4 = False
                            break

                    no_otherK_within_3_4 = True
                    for d in distances_to_otherK:
                        if 0 < d <= 4.0:
                            no_otherK_within_3_4 = False
                            break

                    no_otherCl_within_3_4 = True
                    for d in distances_to_Cl:
                        if d <= 3.4:
                            no_otherCl_within_3_4 = False
                            break

                    if has_SOL_within_3_4 and no_RES_within_3_4 and no_otherK_within_3_4 and no_otherCl_within_3_4:
                        K_with_SOL_within_4_no_RES_no_otherK += 1

                    has_SOL_againwithin_3_4 = False
                    for d in distances_to_SOL:
                        if d <= 3.4:
                            has_SOL_againwithin_3_4 = True
                            break

                    no_RES_within_3_4 = True
                    for d in distances_to_RES:
                        if d <= 3.0:
                            no_RES_within_3_4 = False
                            break

                    no_otherK_within_3_4 = True
                    for d in distances_to_otherK:
                        if 0 < d <= 4.0:
                            no_otherK_within_3_4 = False
                            break

                    has_otherCl_within_3_4 = False
                    for d in distances_to_Cl:
                        if d <= 3.4:
                            has_otherCl_within_3_4 = True
                            break

                    if has_SOL_againwithin_3_4 and no_RES_within_3_4 and no_otherK_within_3_4 and has_otherCl_within_3_4:
                        K_with_SOL_within_4_no_RES_has_otherK_Cl += 1

                    # Count K_with_K_or_Cl_no_SOL_or_RES
                    no_SOL_within_4 = True
                    for d in distances_to_SOL:
                        if d <= 3.7:
                            no_SOL_within_4 = False
                            break

                    no_RES_within_4 = True
                    for d in distances_to_RES:
                        if d <= 3.0:
                            no_RES_within_4 = False
                            break

                    has_otherK_within_3_2 = False
                    for d in distances_to_otherK:
                        if 0 < d <= 3.4:
                            has_otherK_within_3_2 = True
                            break

                    has_Cl_within_3_2 = False
                    for d in distances_to_Cl:
                        if d <= 3.3:
                            has_Cl_within_3_2 = True
                            break

                    if no_SOL_within_4 and no_RES_within_4 and (has_otherK_within_3_2 or has_Cl_within_3_2):
                        K_with_K_or_Cl_no_SOL_or_RES += 1

            result_line = (f"{ts.frame}\t{RES_with_SOL_within_4_no_K_no_RES}\t"
                           f"{RES_with_RES_or_Cl_yes_SOL_or_K}\t{RES_with_RES_or_Cl_no_SOL_or_K}\t{K_with_SOL_within_4_no_RES_no_otherK}\t"
                           f"{K_with_SOL_within_4_no_RES_has_otherK_Cl}\t{K_with_K_or_Cl_no_SOL_or_RES}\n")
            f.write(result_line)
            print(result_line.strip())


if __name__ == "__main__":
    tpr_file = '10.tpr'
    xtc_file = '10.xtc'
    u = mda.Universe(tpr_file, xtc_file)
    analyze_residues(u)
