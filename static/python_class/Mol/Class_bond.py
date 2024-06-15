# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 15:02:54 2023

@author: Toma Legrand
"""

class Bond:
    
    def __init__(self,line_in_file):
        self.line_in_file = line_in_file
        self.atom_num = int((self.line_in_file[3].split())[0])
        self.bond_num = int((self.line_in_file[3].split())[1])
        self.bonds = []
        self.nature_bonds = []
        
    def add_bonds(self):
        """add_bonds: self -> bool
        Reads IDs of two bond-forming atoms from each line of a list 
        of file lines, starting on the line after the atom block. Appends each
        atom pair as a tuple to self.bonds. Returns True, if the bond block 
        length corresponds to the bond number in the header, and False, if the
        bond block is longer or shorter.
        """
        length_bb = True
        start = 4 + self.atom_num
        end = start + self.bond_num
        l = len(self.line_in_file)
        if l - start - 2 == self.bond_num:
            for e in range(start,end):
                line = self.line_in_file[e].split()
                bond = (line[0],line[1])
                self.bonds.append(bond)
        else:
            length_bb = False
        return length_bb
    
    def get_bonds(self):
        """get_bonds: self -> list
        """
        return self.bonds
    
    def build_nature_bonds(self):
        """ build_nature_bonds: self -> list
        Fill the self.nature_bonds list based on the type of bonds between atoms.
        The data from the nature of the bonds is the third item of each line in the bond block.
        """
        start = 4 + self.atom_num
        end = start + self.bond_num
        l = len(self.line_in_file)
        if l - start - 2 == self.bond_num:
            for k in range(start,end):
                line = self.line_in_file[k].split()
                self.nature_bonds.append(line[2])
        return self.nature_bonds


