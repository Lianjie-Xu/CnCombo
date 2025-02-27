# ===============================================================================
#     This file is part of Jwalk.
#     
#     Jwalk - A tool to calculate the solvent accessible surface distance (SASD) 
#     between crosslinked residues.
#     
#     Copyright 2016 Jwalk Inventor and Birkbeck College University of London.
#                          The Jwalk Inventor is: Josh Bullock
# 
# 
#     Jwalk is available under Public Licence.
#     This software is made available under GPL V3
#
#     Please cite your use of Jwalk in published work:
#     
#     J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#     The importance of non-accessible crosslinks and solvent accessible surface distance
#     in modelling proteins with restraints from crosslinking mass spectrometry. 
#     In revision.
#
# ===============================================================================

from numpy import array, zeros, real, sqrt, exp
from collections import OrderedDict
import sys


class Map:
    """ 
    
    A class representing all information from a density map file. 
    NOTE: Currently it can only read the CCP4/MRC  format.
    
    """

    def __init__(self, fullMap, origin, apix, filename, header=[]):
        """
        
        Read a map and its parameters in to Map class instance.
        
        *filename*
            name of map file.
        *origin*    
            origin co-ordinates of the map (x_origin, y_origin, z_origin).
        *apix*
            grid spacing of map.
        *filename*
            filename of the Map instance
            
            NOTE: The *filename* 'build++copy' is reserved for copying of other Map class instances."""
        self.header = header
        self.origin = origin
        self.apix = apix
        self.filename = filename
        self.fullMap = fullMap

    def copy(self):
        """
        
        Return:
            copy of the Map.
        
        """
        copy = Map(self.fullMap.copy(), self.origin[:], self.apix, self.filename, self.header[:])
        return copy

    def box_size(self):
        """
        
        Return:
            size of the map array, in ZYX format.
        
        """
        return self.fullMap.shape

    def x_size(self):
        """
        
        Return:
            x size of the map array in x direction.
        
        """
        return self.fullMap.shape[2]

    def y_size(self):
        """
        
        Return:
            y size of the map array in y direction.
        
        """
        return self.fullMap.shape[1]

    def z_size(self):
        """
        
        Return:
            z size of the map array in z direction.
        
        """
        return self.fullMap.shape[0]


def makeGrid(struct, apix, resolution, filename="None"):
    # Build empty template map based on the size of the protein and the resolution.
    extr = struct.get_extreme_values()
    edge = int(2 * resolution / apix) + 2
    x_size = int((extr[1] - extr[0]) / apix) + edge
    y_size = int((extr[3] - extr[2]) / apix) + edge
    z_size = int((extr[5] - extr[4]) / apix) + edge

    # Origin calculated such that the centre of the map is the centre of mass of the protein.
    x_origin = (extr[1] + extr[0]) / 2 - (apix * x_size / 2.0)
    y_origin = (extr[3] + extr[2]) / 2 - (apix * y_size / 2.0)
    z_origin = (extr[5] + extr[4]) / 2 - (apix * z_size / 2.0)

    newMap = zeros((z_size, y_size, x_size))
    fullMap = Map(newMap, [x_origin, y_origin, z_origin], apix, filename)
    return fullMap


def mapGridPosition(densMap, atom):
    """

    Returns the index of the nearest pixel to an atom, and atom mass (4 values in list form).

    Arguments:

       *densMap*
           Map instance the atom is to be placed on.
       *atom*
           Atom instance.

       """
    origin = densMap.origin
    apix = densMap.apix
    box_size = densMap.box_size()
    x_pos = int(round((atom.x - origin[0]) / apix, 0))
    y_pos = int(round((atom.y - origin[1]) / apix, 0))
    z_pos = int(round((atom.z - origin[2]) / apix, 0))

    if ((densMap.x_size() > x_pos >= 0) and (densMap.y_size() > y_pos >= 0) and (densMap.z_size() > z_pos >= 0)):
        return (x_pos, y_pos, z_pos, atom.mass)
    else:
        return 0


def mark_CA_2types(densMap, prot, aa1, aa2):
    lys_CA = OrderedDict()
    acid_CA = OrderedDict()

    for atom in prot.atomList:

        if atom.res == aa1:

            if atom.atom_name == 'CA':
                pos = mapGridPosition(densMap, atom)
                lys_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

        if atom.res == aa2:
            if atom.atom_name == 'CA':
                pos = mapGridPosition(densMap, atom)
                acid_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

    return lys_CA, acid_CA


def mark_CA_lys_acid(densMap, prot):
    lys_CA = OrderedDict()
    acid_CA = OrderedDict()

    for atom in prot.atomList:

        if atom.res == 'LYS':

            if atom.atom_name == 'CA':
                pos = mapGridPosition(densMap, atom)
                lys_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

        if atom.res_no == 1:
            if atom.atom_name == 'CA':
                pos = mapGridPosition(densMap, atom)
                lys_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

        if atom.res == 'GLU':
            if atom.atom_name == 'CA':
                pos = mapGridPosition(densMap, atom)
                acid_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

        elif atom.res == 'ASP':
            if atom.atom_name == 'CA':
                pos = mapGridPosition(densMap, atom)
                acid_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

    return lys_CA, acid_CA


def mark_CA_UV(densMap, prot, uv_xl):
    xl_list = []
    xl_pair = OrderedDict()
    xl_CA = OrderedDict()
    xl_hack = {}
    xl_pair2 = {}
    c = 0

    with open(uv_xl) as xl_in:
        for line in xl_in:
            col = line.split("|")
            chain1 = col[1].rstrip()
            chain1 = chain1.lstrip()
            chain2 = col[3].rstrip()
            chain2 = chain2.lstrip()
            if chain1 == "x":
                chain1 = " "
            if chain2 == "x":
                chain2 = " "
            xl_list.append([int(col[0]), chain1])
            xl_list.append([int(col[2]), chain2])
            if (int(col[0]), chain1, c) in xl_pair:
                c += 1
            xl_pair[int(col[0]), chain1, c] = [int(col[2]), chain2, c]

    xl_check = []

    for atom in prot.atomList:  # xl order not maintined ...
        if [atom.res_no, atom.chain] in xl_list:
            if atom.atom_name == "CA":
                pos = mapGridPosition(densMap, atom)
                xl_hack[atom.res_no, atom.chain] = (atom.res_no, atom.chain, atom.res)
                xl_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]
                xl_check.append([atom.res_no, atom.chain])

    rem_x = []
    for x in xl_list:
        if x not in xl_check:
            print("ERROR ! Residue", x[0], "-", x[1], "not in pdb structure - please check input files")
            rem_x.append(x)

    for x in xl_pair:
        if [xl_pair[x][0], xl_pair[x][1]] in rem_x:
            del xl_pair[x]
        elif [x[0], x[1]] in rem_x:
            del xl_pair[x]

    for x in xl_pair:
        if (x[0], x[1]) in xl_hack:
            xl_pair2[xl_hack[x[0], x[1]][0], xl_hack[x[0], x[1]][1], xl_hack[x[0], x[1]][2], x[2]] = (xl_pair[x])

    # this is a mess, i'm sorry if you're reading this ... does work though.

    for x in xl_pair2:
        if (xl_pair2[x][0], xl_pair2[x][1]) in xl_hack:
            xl_pair2[x] = [xl_hack[xl_pair2[x][0], xl_pair2[x][1]][0], xl_hack[xl_pair2[x][0], xl_pair2[x][1]][1],
                           xl_hack[xl_pair2[x][0], xl_pair2[x][1]][2], x[3]]

    return xl_pair2, xl_CA


def mark_CA(densMap, prot, lysines=True, acidic=False, cysteines=False):
    lys_CA = OrderedDict()
    acid_CA = OrderedDict()
    cys_CA = OrderedDict()

    for atom in prot.atomList:
        if lysines == True:

            if atom.res == 'LYS':

                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    lys_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

            elif atom.res_no == 1:
                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    lys_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

        elif acidic == True:
            if atom.res == 'GLU':
                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    acid_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

            elif atom.res == 'ASP':
                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    acid_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

        elif cysteines == True:
            if atom.res == 'CYS':
                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    cys_CA[atom.res_no, atom.chain, atom.res] = [pos[0], pos[1], pos[2]]

    return lys_CA, acid_CA, cys_CA


def generate_SAS_2type(densMap, prot, lys_CA, acid_CA):
    inside = {}
    sphere = {}
    radius = {}

    C = 0.8

    radius['CA'] = 1.73 + C
    radius['S'] = 1.67 + C
    radius['N'] = 1.43 + C
    radius['OH'] = 1.30 + C

    for r in radius:

        sphere[r] = []
        rad = int(round(radius[r] / densMap.apix))

        for x in range(-rad, rad + 1):
            for y in range(-rad, rad + 1):
                for z in range(-rad, rad + 1):
                    if (x ** 2 + y ** 2 + z ** 2) <= (rad ** 2):
                        sphere[r].append([x, y, z])

    backbone = ['N', 'CA', 'C', 'O']

    for atom in prot.atomList:
        pos = mapGridPosition(densMap, atom)

        if pos:

            if ((atom.res_no, atom.chain, atom.res) in lys_CA and atom.atom_name not in backbone) or (
                    (atom.res_no, atom.chain, atom.res) in acid_CA and atom.atom_name not in backbone):
                pass
            else:
                if atom.atom_name[:1] == 'C':
                    for (x, y, z) in sphere['CA']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'O':
                    for (x, y, z) in sphere['OH']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'N':
                    for (x, y, z) in sphere['N']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'S':
                    for (x, y, z) in sphere['S']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1

    return densMap


def generate_SAS_UV(densMap, prot, uv):
    inside = {}
    sphere = {}
    radius = {}

    C = 0.8

    radius['CA'] = 1.73 + C
    radius['S'] = 1.67 + C
    radius['N'] = 1.43 + C
    radius['OH'] = 1.30 + C

    for r in radius:

        sphere[r] = []
        rad = int(round(radius[r] / densMap.apix))

        for x in range(-rad, rad + 1):
            for y in range(-rad, rad + 1):
                for z in range(-rad, rad + 1):
                    if (x ** 2 + y ** 2 + z ** 2) <= (rad ** 2):
                        sphere[r].append([x, y, z])

    backbone = ['N', 'CA', 'C', 'O']

    for atom in prot.atomList:
        pos = mapGridPosition(densMap, atom)

        if pos:

            if (atom.res_no, atom.chain, atom.res) in uv and atom.atom_name not in backbone:
                pass
            else:
                if atom.atom_name[:1] == 'C':
                    for (x, y, z) in sphere['CA']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'O':
                    for (x, y, z) in sphere['OH']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'N':
                    for (x, y, z) in sphere['N']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'S':
                    for (x, y, z) in sphere['S']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1

    return densMap


def generate_SAS(densMap, prot, lysines=True, acidic=False, cysteines=False):
    """
    spits off spheres ov radius 2(for the orig atom) + 0.5 * 2(for the CH2 of the cross-linker)
                                        0.5 is to allow for low resolution error.
    slight update 12/4/15:
    CA = 1.73 + 1.68 = 4
    CH2 = 1.68 + 1.68 = 4
    S = 1.67 + 1.68 = 4
    N = 1.43 + 1.68 = 4
    OH = 1.30 + 1.68 = 3
    """

    inside = {}
    sphere = {}
    radius = {}

    C = 0.8

    radius['CA'] = 1.73 + C
    radius['S'] = 1.67 + C
    radius['N'] = 1.43 + C
    radius['OH'] = 1.30 + C

    for r in radius:

        sphere[r] = []
        rad = int(round(radius[r] / densMap.apix))

        for x in range(-rad, rad + 1):
            for y in range(-rad, rad + 1):
                for z in range(-rad, rad + 1):
                    if (x ** 2 + y ** 2 + z ** 2) <= (rad ** 2):
                        sphere[r].append([x, y, z])

    backbone = ['N', 'CA', 'C', 'O']

    for atom in prot.atomList:
        pos = mapGridPosition(densMap, atom)

        if pos:

            if lysines == True and atom.res == 'LYS' and atom.atom_name not in backbone:
                pass
            elif acidic == True and (atom.res == 'GLU' or atom.res == 'ASP') and atom.atom_name not in backbone:
                pass
            elif cysteines == True and atom.res == 'CYS' and atom.atom_name not in backbone:
                pass
            elif lysines == True and atom.res_no == 1 and atom.atom_name not in backbone:
                pass
            else:
                if atom.atom_name[:1] == 'C':
                    for (x, y, z) in sphere['CA']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'O':
                    for (x, y, z) in sphere['OH']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'N':
                    for (x, y, z) in sphere['N']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1
                elif atom.atom_name[:1] == 'S':
                    for (x, y, z) in sphere['S']:
                        if ((densMap.x_size() > (pos[0] + x) >= 0) and (densMap.y_size() > (pos[1] + y) >= 0) and (
                                densMap.z_size() > (pos[2] + z) >= 0)):
                            densMap.fullMap[pos[2] + z][pos[1] + y][pos[0] + x] += 1

    return densMap


def find_empty_space(k, sphere, densMap, CA):
    starters = []
    (x, y, z) = CA[k]
    for (x_s, y_s, z_s) in sphere:
        if ((densMap.x_size() > (x + x_s) >= 0) and (densMap.y_size() > (y + y_s) >= 0) and (
                densMap.z_size() > (z + z_s) >= 0)):
            if densMap.fullMap[z_s + z][y_s + y][x_s + x] <= 0:
                starters.append([x_s + x, y_s + y, z_s + z])
    return starters


def find_surface_voxels(lys_CA, acid_CA, cys_CA, densMap, lysines=True, acidic=False, cysteines=False):
    sphere1 = []
    sphere2 = []
    sphere3 = []
    sphere4 = []
    sphere5 = []

    C = 1.68

    radius = 1.73 + C

    radius4 = int(round(radius / densMap.apix)) + 1  # radius rounded up = 4 with apix = 1
    radius3 = radius4 - 1
    radius2 = radius3 - 1
    radius1 = radius2 - 1
    radius5 = radius4 + 1

    for x in range(-radius1, radius1 + 1):
        for y in range(-radius1, radius1 + 1):
            for z in range(-radius1, radius1 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius1 ** 2):
                    sphere1.append([x, y, z])

    for x in range(-radius2, radius2 + 1):
        for y in range(-radius2, radius2 + 1):
            for z in range(-radius2, radius2 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius2 ** 2):
                    if ([x, y, z]) not in sphere1:
                        sphere2.append([x, y, z])

    for x in range(-radius3, radius3 + 1):
        for y in range(-radius3, radius3 + 1):
            for z in range(-radius3, radius3 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius3 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2:
                        sphere3.append([x, y, z])

    for x in range(-radius4, radius4 + 1):
        for y in range(-radius4, radius4 + 1):
            for z in range(-radius4, radius4 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius4 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3:
                        sphere4.append([x, y, z])

    for x in range(-radius5, radius5 + 1):
        for y in range(-radius5, radius5 + 1):
            for z in range(-radius5, radius5 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius5 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3 and (
                    [x, y, z]) not in sphere4:
                        sphere5.append([x, y, z])

    lys_starters = OrderedDict()
    acid_starters = OrderedDict()
    cys_starters = OrderedDict()

    if lysines == True:
        k_count = 0
        k_buried = 0
        for k in lys_CA:
            k_count += 1
            lys_starters[k] = find_empty_space(k, sphere1, densMap, lys_CA)
            if lys_starters[k] == []:
                lys_starters[k] = find_empty_space(k, sphere2, densMap, lys_CA)
            if lys_starters[k] == []:
                lys_starters[k] = find_empty_space(k, sphere3, densMap, lys_CA)
            if lys_starters[k] == []:
                lys_starters[k] = find_empty_space(k, sphere4, densMap, lys_CA)
            if lys_starters[k] == []:
                lys_starters[k] = find_empty_space(k, sphere5, densMap, lys_CA)
            if lys_starters[k] == []:
                del lys_starters[k]
                k_buried += 1

    if acidic == True:
        k_count = 0
        k_buried = 0
        for k in acid_CA:
            k_count += 1
            acid_starters[k] = find_empty_space(k, sphere1, densMap, acid_CA)
            if acid_starters[k] == []:
                acid_starters[k] = find_empty_space(k, sphere2, densMap, acid_CA)
            if acid_starters[k] == []:
                acid_starters[k] = find_empty_space(k, sphere3, densMap, acid_CA)
            if acid_starters[k] == []:
                acid_starters[k] = find_empty_space(k, sphere4, densMap, acid_CA)
            if acid_starters[k] == []:
                acid_starters[k] = find_empty_space(k, sphere5, densMap, acid_CA)
            if acid_starters[k] == []:
                del acid_starters[k]
                k_buried += 1

    if cysteines == True:
        k_count = 0
        k_buried = 0
        for k in cys_CA:
            k_count += 1
            cys_starters[k] = find_empty_space(k, sphere1, densMap, cys_CA)
            if cys_starters[k] == []:
                cys_starters[k] = find_empty_space(k, sphere2, densMap, cys_CA)
            if cys_starters[k] == []:
                cys_starters[k] = find_empty_space(k, sphere3, densMap, cys_CA)
            if cys_starters[k] == []:
                cys_starters[k] = find_empty_space(k, sphere4, densMap, cys_CA)
            if cys_starters[k] == []:
                cys_starters[k] = find_empty_space(k, sphere5, densMap, acid_CA)
            if cys_starters[k] == []:
                del cys_starters[k]
                k_buried += 1

    return lys_starters, acid_starters, cys_starters


def find_surface_voxels_surface_false(lys_CA, densMap, uv_list=False):
    sphere1 = []
    sphere2 = []
    sphere3 = []
    sphere4 = []

    C = 1.68

    radius = 1.73 + C

    radius4 = int(round(radius / densMap.apix)) + 1  # radius rounded up = 4 with apix = 1
    radius3 = radius4 - 1
    radius2 = radius3 - 1
    radius1 = radius2 - 1

    for x in range(-radius1, radius1 + 1):
        for y in range(-radius1, radius1 + 1):
            for z in range(-radius1, radius1 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius1 ** 2):
                    sphere1.append([x, y, z])

    for x in range(-radius2, radius2 + 1):
        for y in range(-radius2, radius2 + 1):
            for z in range(-radius2, radius2 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius2 ** 2):
                    if ([x, y, z]) not in sphere1:
                        sphere2.append([x, y, z])

    for x in range(-radius3, radius3 + 1):
        for y in range(-radius3, radius3 + 1):
            for z in range(-radius3, radius3 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius3 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2:
                        sphere3.append([x, y, z])

    for x in range(-radius4, radius4 + 1):
        for y in range(-radius4, radius4 + 1):
            for z in range(-radius4, radius4 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius4 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3:
                        sphere4.append([x, y, z])

    lys_starters = OrderedDict()

    tell_buried = []

    k_count = 0
    k_buried = 0
    for k in lys_CA:
        k_count += 1
        lys_starters[k] = find_empty_space(k, sphere1, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere2, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere3, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere4, densMap, lys_CA)
        if lys_starters[k] == []:
            tell_buried.append(k)
            del lys_starters[k]
            k_buried += 1

    rem_x = []
    if len(tell_buried) > 0 and uv_list == True:
        print("ERROR - ", k_buried, " buried residue(s) in xl_list:")
        for t in tell_buried:
            print(str(t[0]) + "-" + str(t[1]) + "-" + str(t[2]))
            rem_x.append([t[0], t[1]])
        # sys.exit(2)
    return lys_starters, rem_x


def find_surface_voxels2(CA, densMap):
    sphere1 = []
    sphere2 = []
    sphere3 = []
    sphere4 = []
    sphere5 = []
    sphere6 = []

    C = 1.68

    radius = 1.73 + C

    radius4 = int(round(radius / densMap.apix)) + 1  # radius rounded up = 4 with apix = 1
    radius3 = radius4 - 1
    radius2 = radius3 - 1
    radius1 = radius2 - 1
    radius5 = radius4 + 1
    radius6 = radius5 + 1

    for x in range(-radius1, radius1 + 1):
        for y in range(-radius1, radius1 + 1):
            for z in range(-radius1, radius1 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius1 ** 2):
                    sphere1.append([x, y, z])

    for x in range(-radius2, radius2 + 1):
        for y in range(-radius2, radius2 + 1):
            for z in range(-radius2, radius2 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius2 ** 2):
                    if ([x, y, z]) not in sphere1:
                        sphere2.append([x, y, z])

    for x in range(-radius3, radius3 + 1):
        for y in range(-radius3, radius3 + 1):
            for z in range(-radius3, radius3 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius3 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2:
                        sphere3.append([x, y, z])

    for x in range(-radius4, radius4 + 1):
        for y in range(-radius4, radius4 + 1):
            for z in range(-radius4, radius4 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius4 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3:
                        sphere4.append([x, y, z])

    for x in range(-radius5, radius5 + 1):
        for y in range(-radius5, radius5 + 1):
            for z in range(-radius5, radius5 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius5 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3 and (
                    [x, y, z]) not in sphere4:
                        sphere5.append([x, y, z])

    for x in range(-radius6, radius6 + 1):
        for y in range(-radius6, radius6 + 1):
            for z in range(-radius6, radius6 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius6 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3 and (
                    [x, y, z]) not in sphere4 and ([x, y, z]) not in sphere5:
                        sphere6.append([x, y, z])

    starters = OrderedDict()

    k_count = 0
    k_buried = 0
    for k in CA:
        k_count += 1
        starters[k] = find_empty_space(k, sphere1, densMap, CA)
        if starters[k] == []:
            starters[k] = find_empty_space(k, sphere2, densMap, CA)
        if starters[k] == []:
            starters[k] = find_empty_space(k, sphere3, densMap, CA)
        if starters[k] == []:
            starters[k] = find_empty_space(k, sphere4, densMap, CA)
        if starters[k] == []:
            starters[k] = find_empty_space(k, sphere5, densMap, CA)
        if starters[k] == []:
            del starters[k]
            k_buried += 1

    return starters


def find_surface_voxels_2type(lys_CA, acid_CA, densMap):
    sphere1 = []
    sphere2 = []
    sphere3 = []
    sphere4 = []
    sphere5 = []
    sphere6 = []

    C = 1.68

    radius = 1.73 + C

    radius4 = int(round(radius / densMap.apix)) + 1  # radius rounded up = 4 with apix = 1
    radius3 = radius4 - 1
    radius2 = radius3 - 1
    radius1 = radius2 - 1
    radius5 = radius4 + 1
    radius6 = radius5 + 1

    for x in range(-radius1, radius1 + 1):
        for y in range(-radius1, radius1 + 1):
            for z in range(-radius1, radius1 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius1 ** 2):
                    sphere1.append([x, y, z])

    for x in range(-radius2, radius2 + 1):
        for y in range(-radius2, radius2 + 1):
            for z in range(-radius2, radius2 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius2 ** 2):
                    if ([x, y, z]) not in sphere1:
                        sphere2.append([x, y, z])

    for x in range(-radius3, radius3 + 1):
        for y in range(-radius3, radius3 + 1):
            for z in range(-radius3, radius3 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius3 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2:
                        sphere3.append([x, y, z])

    for x in range(-radius4, radius4 + 1):
        for y in range(-radius4, radius4 + 1):
            for z in range(-radius4, radius4 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius4 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3:
                        sphere4.append([x, y, z])

    for x in range(-radius5, radius5 + 1):
        for y in range(-radius5, radius5 + 1):
            for z in range(-radius5, radius5 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius5 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3 and (
                    [x, y, z]) not in sphere4:
                        sphere5.append([x, y, z])

    for x in range(-radius6, radius6 + 1):
        for y in range(-radius6, radius6 + 1):
            for z in range(-radius6, radius6 + 1):
                if (x ** 2 + y ** 2 + z ** 2) <= (radius6 ** 2):
                    if ([x, y, z]) not in sphere1 and ([x, y, z]) not in sphere2 and ([x, y, z]) not in sphere3 and (
                    [x, y, z]) not in sphere4 and ([x, y, z]) not in sphere5:
                        sphere6.append([x, y, z])

    lys_starters = OrderedDict()
    acid_starters = OrderedDict()

    k_count = 0
    k_buried = 0
    for k in lys_CA:
        k_count += 1
        lys_starters[k] = find_empty_space(k, sphere1, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere2, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere3, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere4, densMap, lys_CA)
        if lys_starters[k] == []:
            lys_starters[k] = find_empty_space(k, sphere5, densMap, lys_CA)
        if lys_starters[k] == []:
            del lys_starters[k]
            k_buried += 1

    for k in acid_CA:
        k_count += 1
        acid_starters[k] = find_empty_space(k, sphere1, densMap, acid_CA)
        if acid_starters[k] == []:
            acid_starters[k] = find_empty_space(k, sphere2, densMap, acid_CA)
        if acid_starters[k] == []:
            acid_starters[k] = find_empty_space(k, sphere3, densMap, acid_CA)
        if acid_starters[k] == []:
            acid_starters[k] = find_empty_space(k, sphere4, densMap, acid_CA)
        if acid_starters[k] == []:
            acid_starters[k] = find_empty_space(k, sphere5, densMap, acid_CA)
        if acid_starters[k] == []:
            del acid_starters[k]
            k_buried += 1

    return lys_starters, acid_starters


def remove_duplicates(final_xl, uv_list=False):
    if uv_list == True:
        return final_xl

    hold_list = []
    hold = {}
    done = []
    filtered = {}
    filter_list = []
    for (k1, k2, sasd) in final_xl:

        if [k1, k2] not in hold_list:
            hold_list.append([k1, k2])

            hold[k1, k2] = sasd


        elif [k1, k2] in hold_list:
            if (sasd < hold[k1, k2]):
                filtered[(k1, k2, sasd)] = final_xl[(k1, k2, sasd)]
                filter_list.append([k1, k2])
            else:
                filtered[(k1, k2, hold[k1, k2])] = final_xl[(k1, k2, hold[k1, k2])]
                filter_list.append([k1, k2])

    for (k1, k2, sasd) in final_xl:
        if [k1, k2] not in filter_list:
            filter_list.append([k1, k2])
            filtered[(k1, k2, sasd)] = final_xl[(k1, k2, sasd)]

    return filtered


def remove_crosslinks(xl_pair, rem_x):
    delly = []
    for x in xl_pair:
        if [xl_pair[x][0], xl_pair[x][1]] in rem_x:
            delly.append(x)
        elif [x[0], x[1]] in rem_x:
            delly.append(x)
    for d in delly:
        del xl_pair[d]

    return xl_pair
