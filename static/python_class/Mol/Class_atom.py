# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 19:03:31 2023

@author: Toma Legrand
"""

class Atom:
    
    def __init__(self,line_in_file):
        self.line_in_file = line_in_file
        self.atom_num = int((self.line_in_file[3].split())[0])
        self.atoms = {}
        
    def add_atoms(self):
        """add_atoms: self
        Reads atom IDs, atom elements and atom xyz coordinates from a list 
        of file lines, starting on line 5. Adds atoms to self.atoms, using the 
        IDs as keys and [elements, (x,y,z)] as values of each key.
        """
        for num in range(4,4+self.atom_num):
            index = num-3
            line = self.line_in_file[num].split()
            atom_type = line[3]
            atom_xyz = (float(line[0]),float(line[1]),float(line[2]))
            self.atoms[index] = []
            self.atoms[index].append(atom_type)
            self.atoms[index].append(atom_xyz)
    
    def get_atoms(self):
        """get_atoms: self -> dict
        """
        return self.atoms

            
    

    
    
    
            
        

    