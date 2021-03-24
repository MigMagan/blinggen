#! /usr/bin/env python

import textwrap

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class cell:
    """
    This class is a MCNP cell constructed with its surfaces.

    A surface is identified as a number, and the cells are constructed 
    recursively using AND (blankspace) OR (:) logic. So far, this is 
    just to be able to write cells 
    """
    
    def __init__ (self, surfaces, n=0, logic = "and"):
        self.surfaces = [surf for surf in surfaces]  #Array describing the surfaces
        self.n = n  # Cell number
        if logic == "and":
            self.logic = ["and" for i in surfaces[1:]]
        elif logic == "or":
            self.logic = ["or" for i in surfaces[1:]]
        else:
            print("logic not understood")
            pass
    
    def __add__ (self, o):
        result = cell([0])
        if all(a == "or" for a in self.logic):
            result.logic = self.logic +["or"]
            result.surfaces = self.surfaces + [o]
        
        else:
            result.logic = ["or"]
            result.surfaces = [self, o]
            result.n = self.n
        return result

    def __neg__ (self):
        result = cell([0], -self.n)
        result.surfaces = [-surf for surf in self.surfaces]
        result.logic = ["and" if i == "or" else "or" for i in self.logic]
        return result
    
    def __sub__(self, o):
        s = -o
        return self * s
    
    def __mul__(self, o):  # Intersection
        result = cell([0])
        if all(a == "and" for a in self.logic):
            result.logic = self.logic + ["and"]
            result.surfaces = self.surfaces + [o]

        else:
            result.logic = ["and"]
            result.surfaces = [self, o]
        result.n = self.n
        return result 
# ========================== END OF CLASS DEFINITION ========================== 

def __surfcard(cell):
    "Internal method to build the surface part of the mcnp card"
    card = "("
    for i, surf in enumerate(cell.surfaces):
        if type(surf) == int:
            card+=(str(surf)) 
        else:
            card+=__surfcard(surf)
        if i != len(cell.surfaces)-1:
            if cell.logic[i] == "and":
                card+=" "
            elif cell.logic[i] == "or":
                card+=":"
    card+=")" 
    return card
   
def mcnpcard(cell, mat=0, ro=None):
    " dump the mcnp surface definition of the cell"
    uwcard = "{0}  {1}  {2} ".format(cell.n, mat, str(ro or ''))
    uwcard = uwcard+__surfcard(cell)+"\n      IMP:N=1"
    card = "\n".join(textwrap.wrap(uwcard, width=75, subsequent_indent='     '))
    return card

