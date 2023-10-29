# Importing libraries
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition

# Functions
def get_nn_nos(nn_list):
    nn_nos = []
    for nn in nn_list:
        nn_nos.append(len(nn))
    return nn_nos

def get_indices(no_list, nn_nos):
    """
    This function returns indices of sites with specified number of nearest neighbour.
    
    no_list : a list containing the number of nearest neighbor to consider.
    nn_nos : nearest neibour number list
    """
    ind = []
    for n in no_list:
        for i in range(len(nn_nos)):
            if n == nn_nos[i]:
                ind.append(i)
    return ind

def center_of_mass(poscar):
    """
    Returns the center of mass of a structure
    """
    center = np.zeros(3)
    total_weight = 0
    for site in poscar.structure.sites:
        wt = site.species.weight
        center += site.coords * wt
        total_weight += wt
        
    return center / total_weight

def find_point_d_from_pos2(pos1, pos2, d=1):
    """
    Gives the a for the quadratic equation at^2 + bt + c = 0
    Equation of the line connecting pos1 (center of mass) and pos2 is
    
    x3 = (x2-x1) + x1
    y3 = (y2-y1) + y1
    z3 = (z2-z1) + z1
    
    We want to determine a points on the line with d distance from the second point.
    
    d is the distance
    """
    # Calculate a
    a = (pos2[0]-pos1[0])**2 + (pos2[1]-pos1[1])**2 + (pos2[2]-pos1[2])**2
    t_in = 1 - d/np.sqrt(a)
    t_out = 1 + d/np.sqrt(a)
    
    in_pos = [(pos2[0]-pos1[0])*t_in + pos1[0], (pos2[1]-pos1[1])*t_in + pos1[1], (pos2[2]-pos1[2])*t_in + pos1[2]]
    out_pos = [(pos2[0]-pos1[0])*t_out + pos1[0], (pos2[1]-pos1[1])*t_out + pos1[1], (pos2[2]-pos1[2])*t_out + pos1[2]]
    return [in_pos, out_pos]


def add_pass_el(poscar, com, ind, d=1, pass_el="H", filename="./POSCAR2"):
    pos1 = com
    poscar2 = poscar.structure
    for i in ind:
        # print(i)
        pos2 = poscar.structure.sites[i].coords
        in_pos, out_pos = find_point_d_from_pos2(pos1, pos2,d=d)
    # pos2, in_pos, out_pos
        poscar2.append(species=Composition(pass_el), coords=out_pos, coords_are_cartesian=True)
    
    poscar2 = Poscar(poscar2)
    poscar2.write_file(filename)
    return poscar2

def find_distance(pos1, pos2):
    return np.sqrt((pos2[0]-pos1[0])**2 + (pos2[1]-pos1[1])**2 + (pos2[2]-pos1[2])**2)

def passivate_element_to_outer_layer(radius, no_nn_to_consider, pass_el="H", d=1, filename1 = "./POSCAR", filename2 = "./POSCAR", filename3 = "./POSCAR"):
    poscar = Poscar.from_file(filename1) # Reading the original POSCAR of the nanodots
    com = center_of_mass(poscar) # Get the center of mass
    nn_list = poscar.structure.get_all_neighbors(radius) # List of nearest neighbour of all the elements
    nn_nos = get_nn_nos(nn_list) # Calculating the number of nearest neighbour for each atoms
    
    print("Radius considered for finding nearest neighbour (in Angstrom) : ", radius)
    print("Minimum no of nearest neighbours : ", min(nn_nos))
    print("Maximum no of nearest neighbours : ", max(nn_nos))
    print("Atoms with following number of are considered outer atoms : ", no_nn_to_consider)
    
    ind = get_indices(no_nn_to_consider, nn_nos) # Get the indices of the outer atom

    # Creating a POSCAR with just the outer atoms
    lattice = poscar.structure.lattice
    species = [poscar.structure.sites[i].species for i in ind]
    coords = [poscar.structure.sites[i].coords for i in ind]

    poscar1 = Structure(lattice=lattice, species=species, coords=coords, coords_are_cartesian=True)
    poscar1 = Poscar(poscar1)
    poscar1.write_file(filename2) # POSCAR containing the outer atoms

    # Poscar with hydrogen added to the sides of the nanostructure
    poscar2 = add_pass_el(poscar, com, ind, d=d, pass_el=pass_el, filename=filename3)
    
    return poscar2
    
if __name__ == "__main__":
    # Inputs from the user
    radius = 3 # Radius around an atom for calculating the nearest neighbour
    no_nn_to_consider = [1,2,3] # Atoms with number of atoms to be considered as outer atoms
    pass_el = "H" # Passivating element
    d = 1 # Hydrogen to be added at a distance of 1 Angstrom from the outer atoms

    # Filenames
    filename1 = "./POSCAR" # To read the Original POSCAR of the nanodots
    filename2 = "./POSCAR1" # To write the POSCAR containing the outer atoms
    filename3 = "./POSCAR2" # To write the POSCAR containing nanodots passivated with hydrogen
    poscar2 = passivate_element_to_outer_layer(radius=radius, no_nn_to_consider=no_nn_to_consider, pass_el=pass_el, d=d, filename1=filename1, filename2=filename2, filename3=filename3)
    
