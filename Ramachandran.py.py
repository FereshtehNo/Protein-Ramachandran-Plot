import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from Bio.PDB import PDBParser, PPBuilder
import warnings
from Bio import BiopythonWarning

def calculate_phi_psi(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    ppb = PPBuilder()
    model = structure[0]  # Assuming single model in the PDB file

    phi_values = []
    psi_values = []

    for chain in model:
        polypeptides = ppb.build_peptides(chain)
        for polypeptide in polypeptides:
            phi_psi = polypeptide.get_phi_psi_list()
            for residue_phi_psi in phi_psi:
                phi, psi = residue_phi_psi
                if phi is not None and psi is not None:  # Exclude None values
                    phi_values.append(np.rad2deg(phi))  # Convert phi to degrees
                    psi_values.append(np.rad2deg(psi))  # Convert psi to degrees

    return phi_values, psi_values

def plot_ramachandran(phis, psis):
    plt.figure(figsize=(9, 9))
    ax = plt.gca()

    # Define regions
    allowed_region = ((-180, -30), (-180, 180))
    favored_region = ((-180, 30), (-180, 180))

    # Set background colors
    ax.set_facecolor('#FFF6F6')  # Mild pink
    plt.plot((-180, 180), (-180, 180), color='#FEFDD4')  # Mild yellow

    # Convert psi and phi angles to radians
    phis_rad = np.radians(phis)
    psis_rad = np.radians(psis)

    # Plot points
    plt.scatter(phis_rad, psis_rad, c=phis_rad, cmap='rainbow', s=5, alpha=0.7)

    # Set plot limits and labels
    plt.xlim(-np.pi, np.pi)
    plt.ylim(-np.pi, np.pi)
    plt.xticks(np.radians(range(-180, 181, 45)), range(-180, 181, 45))
    plt.yticks(np.radians(range(-180, 181, 45)), range(-180, 181, 45))
    plt.xlabel('Phi')
    plt.ylabel('Psi')
    plt.title('Ramachandran Plot')

    # Create colorbar
    cbar = plt.colorbar(orientation='vertical', aspect=40, pad=0.04)
    cbar.set_label('Degrees')

    # Add regions
    for phi, psi in zip(phis, psis):
        if (phi >= allowed_region[0][0] and phi <= allowed_region[0][1] and
                psi >= allowed_region[1][0] and psi <= allowed_region[1][1]):
            color = 'green'  # Allowed region
            radius = 0.5
        elif (phi >= favored_region[0][0] and phi <= favored_region[0][1] and
              psi >= favored_region[1][0] and psi <= favored_region[1][1]):
            color = 'orange'  # Favored region
            radius = 0.2
        else:
            color = 'red'  # Other regions
            radius = 0.1

        # Plot circles around allowed and favored regions
        if color in ['green', 'orange']:
            circle = Circle((phi, psi), radius, color=color, fill=False)
            ax.add_patch(circle)

    # Create legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=5, label='Allowed'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=5, label='Favored'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=5, label='Other')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.grid(True)
    plt.show()

# Suppress PDBConstructionWarning
warnings.simplefilter('ignore', BiopythonWarning)

# Example usage
pdb_file = "D:/Phd-classes/Final-project/ProcessedData/4YDF.pdb"
phis, psis = calculate_phi_psi(pdb_file)
plot_ramachandran(phis, psis)