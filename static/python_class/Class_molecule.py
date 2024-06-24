# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 19:32:24 2023

@author: Toma Legrand
"""

from static.python_class.Mol.Class_bond import Bond
from static.python_class.Mol.Class_atom import Atom
from static.python_class.Mol.Class_file import Mol_File
import re
from rdkit import Chem

class Molecule:
    
    def __init__(self):
        self.line_in_file = []
        self.atoms = {}
        self.bonds = {}
        self.molecule_graph = {}
     
    def build_graph(self):
        """build_graph: self -> dict
        Constructs and returns dict self.molecule_graph. Inserts the atom 
        indices from self.atoms as keys and, if the key index is part of a 
        bond pair in self.bonds, assigns the other index from that bond pair
        as a value to the key.
        """
        self.molecule_graph = {}
        for key in self.atoms.keys():
            self.molecule_graph[key] = []
            for pair in self.bonds:
                if str(key) in pair:
                    for i in pair:
                        if i != str(key) and i not in self.molecule_graph[key]:
                            self.molecule_graph[key].append(int(i))
        return self.molecule_graph
        
#For command LOAD_VERIFY
    def load_file(self,file_str):
        """load_file: self x string -> bool
        Read file with file_name. 
        Splits file into a list of lines and assigns it to self.line_in_file.
        """
        self.line_in_file = file_str.splitlines()
        line_one = bool(re.match("[A-Z][a-z]?[0-9]*",str(self.line_in_file[0].split()[0])))
        if line_one == False or len(self.line_in_file) <= 4:
            return False
        print(self.line_in_file)
        
    def initialize_variables(self):
        """initialize_variables: self
        Assigns self.atoms and self.bonds to the dict object.atoms from the 
        class Atom and the list object.bonds from the class Bond.
        """
        atoms = Atom(self.line_in_file)
        atoms.add_atoms()
        self.atoms = atoms.get_atoms()
        bonds = Bond(self.line_in_file)
        bonds.add_bonds()
        self.bonds = bonds.get_bonds()
        
    def verification(self,atom_list):
        """verification: self x list -> int
        Checks, whether the methods .verify_organic(atom_list), .add_bonds()
        and .verify_atoms are True. Returns 4, if all methods return True.
        Returns 1, 2, 3 if one of the methods returns False.
        """
        file = Mol_File(self.atoms,self.bonds)
        bond = Bond(self.line_in_file)
        if file.verify_organic(atom_list):
            if bond.add_bonds():
                if file.verify_atoms():
                    return 4    #verification successful
                else:
                    return 3    #atoms not entirely covered by bond block
            else:
                return 2        #bond block length doesn't correspond to bond num
        else:
            return 1            #atom not in atom_list (not organic)
               
        
#For command COUNT_ELEMENTS     
    def count(self,element):
        """count: self x str -> int
        Returns number of element in the dict self.atoms.
        """
        number = 0
        for key in self.atoms.keys():
            if self.atoms[key][0] == element:
                number += 1
        return number

#For command 3D_DISTANCE    
    def coordinates_sum(self):
        """coordinate_sum: self -> int
        Adds 1 to coord_not_zero for each xyz coordinate in self.atoms unequal
        0.0. If all coordinates are 0.0, returns 0.
        """
        coord_not_zero = 0
        for key in self.atoms.keys():
            for coord in self.atoms[key][1]:
                if coord != 0.0:
                    coord_not_zero += 1
        return coord_not_zero
    
    def room_distance(self,atom_one,atom_two):
        """room_distance: self x int x int -> string or float
        Determines the 3D distance between two atom indices, rounded to two 
        positions after the comma, from the xyz coordinates in self.atoms[atom_one]
        and self.atoms[atom_two]. Returns error string, if one atom index is
        not present in self.atoms.keys() or self.atoms[key] does not contain
        any coordinates.
        """
        import math
        if atom_one not in self.atoms.keys():
            return "Wrong ID"
        elif atom_two not in self.atoms.keys():
            return "Wrong ID"
        elif self.atoms[atom_one][1] == ():
            return "No XYZ"
        elif self.atoms[atom_two][1] == ():
            return "No XYZ"
        else:
            o = self.atoms[int(atom_one)][1]
            t = self.atoms[int(atom_two)][1]
            distance = math.sqrt((float(t[0])-float(o[0]))**2+(float(t[1])-float(o[1]))**2+(float(t[2])-float(o[2]))**2)
            distance = round(distance,2)
            return distance
        
#For command ATOM_NEIGHBOURS                           
    def get_neighbour_num(self,atom):
        """get_neighbour_num: self x index -> int or bool
        If the atom index provided is a key in self.molecule_graph, returns the
        number of values (indices of connected atoms) for that atom index. 
        Returns False, if atom not in self.molecule_graph.keys().
        """
        if atom in self.molecule_graph.keys():
            neighbour_num = len(self.molecule_graph[atom])
            return neighbour_num
        else:
            return False
  
#For command 2D_DISTANCE
    def path_distance(self,atom_one,atom_two):
        """path_distance: self x int x int -> int or str
        If both provided atom indices are keys in self.molecule_graph, returns
        the number of edges (bonds) between the two atom vertices. If atom_one
        and atom_two are identicall, returns 0. If one of the provided atoms
        is not a key in self.molecule_graph, returns "False".
        """
        if atom_one in self.molecule_graph.keys() and atom_two in self.molecule_graph.keys():
            explored = [atom_one]
            waiting_queue = [atom_one]
            distance = {}
            distance[atom_one] = 0
            while atom_two not in waiting_queue:
                v = waiting_queue.pop(0)
                for w in self.molecule_graph[v]:
                    if w not in explored:
                        waiting_queue.append(w)
                        explored.append(w)
                        distance[w] = distance[v] + 1
            return distance[atom_two]
        elif atom_one == atom_two:
            return 0
        else:
            return "False"
  
#For command FIND_RING       
    def find_ring(self,atom_list):
        """find_ring: self x list -> int
        Return the results of the number of insaturation in the molecule minus the number of 
        multiple bonds in it (double/triple bonds). If there is no insaturation in the molecule, 
        return 0 else, return the result of the operation. This result correspond to the number of 
        ring in the molecule.
        """
        count = 0
        bond = Bond(self.line_in_file)
        file = Mol_File(self.atoms,self.bonds)
        multiple_bond = bond.build_nature_bonds()
        for k in multiple_bond:
            if int(k) == 2 or int(k) == 3:
                count += 1
        nb_insaturation = file.calcul_insaturation(atom_list)
        if nb_insaturation - count == 0:
            return 0
        else:
            return int(nb_insaturation-count)

# For command TO_SMILES
    def to_smiles(self):
        mol_block = "\n".join(self.line_in_file)
        mol = Chem.MolFromMolBlock(mol_block)
        if mol:
            return Chem.MolToSmiles(mol)
        else:
            return None


        