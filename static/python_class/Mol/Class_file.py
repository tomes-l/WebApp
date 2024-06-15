# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 19:45:31 2023

@author: Toma Legrand
"""

class Mol_File:
    
    def __init__(self,atom_block,bond_block):
        self.atom_block = atom_block
        self.bond_block = bond_block
    
    def verify_organic(self,atom_list):
        """verify_organic: self x list -> bool
        Returns True, if all atom elements in self.atom_block are in 
        atom_list. Returns False, if at least one atom element in 
        self.atom_block is not in atom_list.
        """
        organic = True
        for key in self.atom_block.keys():
            if self.atom_block[key][0] not in atom_list:
                organic = False
                break
        return organic
    
    def verify_atoms(self):
        """verify_atoms: self -> bool
        Returns True, if all atom IDs in self.atom_block are contained in
        elements of self.bond_block. Returns False, if at least one atom ID in 
        self.atom_block is not in elements of self.bond_block.
        """
        all_atoms = []
        for pair in self.bond_block:
            all_atoms.append(int(pair[0]))
            all_atoms.append(int(pair[1]))
        atom_in_bond = True
        for key in self.atom_block.keys():
            if key not in all_atoms:
                atom_in_bond = False
                break
        return atom_in_bond
    
    def calcul_insaturation(self,atom_list):
        """calcul_insaturation: self x list -> int
        Return the number of insaturation of the molecule based on the type of 
        each atoms of a line in the atom block. If there is a division of 0 in the formula, return 0.
        Else, return the results of the formula.
        """
        count_C = 0
        count_H = 0
        count_N = 0
        count_X = 0
        for key in self.atom_block.keys():
            if self.atom_block[key][0] == atom_list[2]:
                count_C += 1
            elif self.atom_block[key][0] == atom_list[0]:
                count_H += 1
            elif self.atom_block[key][0] == atom_list[3] or self.atom_block[key][0] == atom_list[5]:
                count_N += 1
            elif self.atom_block[key][0] == atom_list[7] or self.atom_block[key][0] == atom_list[8] or self.atom_block[key][0] == atom_list[9] or self.atom_block[key][0] == atom_list[10]:
                count_X += 1
        validation_zero = 2*count_C+2-count_H+count_N-count_X
        if validation_zero == 0:
            return 0
        else:
            return int(validation_zero/2)
            
        
           