from openbabel import openbabel, pybel
from scipy.spatial import distance
import numpy as np

elements = ["Bq","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Te", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm","Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

class MZMol(pybel.Molecule):
    def __init__(self, mol: openbabel.OBMol, parser_func=None):
        super(MZMol, self).__init__(mol)
        self._parser_func = parser_func
        self._coords = np.array([atom.coords for atom in self.atoms])

    """
    Atoms related properties
    """
    @property
    def natom(self):
        """
        Get number of atoms in molecule.
        """
        return len(self.atoms)

    @property
    def Zs(self):
        """
        Get atomic number array.
        """
        return np.array([atom.atomicnum for atom in self.atoms])

    @property
    def symbols(self):
        """
        Get chemical symbol array.
        """
        return np.array([elements[Z] for Z in self.Zs])

    @property
    def masses(self):
        """
        Get mass array. (in amu)
        """
        return np.array([atom.atomicmass for atom in self.atoms])

    """
    Coordinates related properties
    """
    @property
    def coords(self):
        """
        Export molecular coordinates in xyz format. (in Angstrom)
        """
        return self._coords

    @coords.setter
    def coords(self, new_coords):
        self._coords = new_coords

    @property
    def xyz_str(self):
        """
        Export this molecule to string in xyz file format. (in Angstrom)
        """
        return self.write(format="xyz")

    @property
    def distance_matrix(self):
        """
        Calculate molecular coordinates in distance matrix format. (in Angstrom)
        Take H2O molecule as an example:
                    1          2          3
                    H          O          O
        1  H    0.000000   0.962698   1.530008
        2  O    0.962698   0.000000   0.962698
        3  H    1.530008   0.962698   0.000000
        """
        return distance.squareform(self.distance_vector)

    @property
    def distance_vector(self):
        """
        Convert distance matrix to vector. (in Angstrom)
        Take H2O molecule as an example:
            r12       r13       r23
        [0.962698, 1.530008, 0.962698]
        """
        return distance.pdist(self.coords)
    
    @property
    def J_r_xyz(self):
        """
        Calculate the Jacobi matrix
        """
        J_r_xyz = np.zeros((self.distance_vector.shape[0], 3 * self.natom))
        k = 0
        for i in range(0, self.natom - 1):
            for j in range(i + 1, self.natom):
                J_r_xyz[k, (3 * i): (3 * i + 3)] = (self.coords[i] - self.coords[j]) / self.distance_vector[k]
                J_r_xyz[k, (3 * j): (3 * j + 3)] = -1.0 * J_r_xyz[k, (3 * i): (3 * i + 3)]
                k += 1
                
        return J_r_xyz

if __name__ == "__main__":
    from openbabel import pybel
    _mol = pybel.readstring(format="xyz", string="""5
    
C   -0.37080090   -0.13842560   -0.26234460
H    0.52405100   -0.92893340   -0.16449230
H   -0.75608510    0.09429450    0.62230030
H    0.06355260    0.73483450   -0.67571110
H   -1.14601250   -0.54180340   -1.11155510""")

    mol = MZMol(_mol)
    print("=====distance matrix=====")
    print(mol.distance_matrix)
    print("=====distance vector=====")
    print(mol.distance_vector)

    np.savetxt("J_r_xyz.txt", mol.J_r_xyz)
    np.savetxt("J_r_xyz_2.txt", mol.J_r_xyz_2)