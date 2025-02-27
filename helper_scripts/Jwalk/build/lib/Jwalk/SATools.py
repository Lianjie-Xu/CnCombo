#===============================================================================
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
#===============================================================================

import sys
def expand_points(atom,sphere):
    
    points = {}
    atom.point_list = []
    CH2 = 1.68
    radius = {"N":1.43,
              "O":1.30,
              "C":1.68,
              "S":1.67,
              }
    
    r = radius[atom.atom_name[0]] + CH2
    (x,y,z) = (atom.x, atom.y, atom.z)
    
    for s in sphere:
        x1 = x + s[0]*r
        y1 = y + s[1]*r
        z1 = z + s[2]*r
        points[x1,y1,z1] = 0

    return points

def create_unit_sphere():
    # creates 30 points
    
    from math import cos, sin
    unit_sphere = []
    
    unit_sphere.append( [0.0, 1.0, 0.0] )
    unit_sphere.append( [0.0, -1.0, 0.0] )
    
    nstep = 5
    PI = 3.1415926536
    theta = PI/nstep
    arc = theta
        
    for istep in range(nstep):
        istep = istep + 1    # to change range from 0--9 to 1--10
        y1 = cos(istep*theta)
        r2 = sin(istep*theta)
        ndot2= 2*PI*r2/arc   # the circumference at that radius / proportion of pi
        if ndot2 == 0.0:
            continue
        theta2 = 2*PI/ndot2
        for idot in range(int(ndot2)):
            idot = idot + 1  # to change range from 0-- to 1--
            x2 = r2*cos(idot*theta2)
            z2 = r2*sin(idot*theta2)
            unit_sphere.append( [x2, y1, z2] )
    
    return unit_sphere 
	
def calculate_if_SA(prot,ROI,uv = False):
    
    string_dict = {"LYS":"lysines",
                "CYS":"cysteines",
                "ASP":"acidic residues",
                "GLU":"acidic residues", 
                "VAL":"valines",
                "ILE":"isoleucines",
                "LEU":"leucines",
                "ARG":"arginines",
                "PRO":"prolines",
                "GLY":"glycines",
                "ALA":"alanines",
                "TRP":"tryptophans",
                "PHE":"phenylalanines",
                "SER":"serines",
                "GLN":"glutamines",
                "HIS":"histidines",
                "MET":"methionines",
                "THR":"threonines",
                "ASN":"asparagines",
                "TYR":"tyrosines"
                }
    
    radius = {"N":1.43,
              "O":1.30,
              "C":1.68,
              "S":1.67,
              }
    sd_res = False
    sphere = create_unit_sphere()
    
    SA_res = {}
    CH2 = 1.68
    for atom in prot.atomList:
        atom.clash = 1
        if (atom.res_no,atom.chain,atom.res) in ROI:
            sd_res = atom.res
            points = expand_points(atom,sphere)
            
            for p in points:
                (x,y,z) = (p[0],p[1],p[2])
                for atom2 in prot.atomList:
                    if atom2.res != atom.res or (atom2.res == atom.res and atom2.atom_name != atom.atom_name):
                        r = radius[atom2.atom_name[0]] + CH2
                        #need to transpose x,y,z
                        (tx,ty,tz) = (x-atom2.x,y-atom2.y,z-atom2.z)
                        if tx**2 + ty**2 + tz**2 <= r**2:
                            points[p] = 1
                            break
                if points[p] == 0:
                    atom.clash = 0
                    break
        
            if atom.clash == 0:
                SA_res[atom.res_no,atom.chain,atom.res] = ROI[atom.res_no,atom.chain,atom.res]
    
    if sd_res == False:
        print ("ERROR: Please type amino acid in three letter code format")
        sys.exit(2)
    
    if uv == True:
        if len(ROI) != len(SA_res):
            print ("the following residues are not solvent accessible")
            for s in SA_res:
                if s not in ROI:
                    print ("%s-%s-%s" % (s[2],s[0],s[1]))
    
    elif sd_res == "LYS":
        print ("%d %s and 1 N-terminus of which %d are on the surface" % (len(ROI)-1,string_dict[sd_res], len(SA_res)))
                    
    else:            
        print ("%d %s of which %d are on the surface" % (len(ROI),string_dict[sd_res], len(SA_res)))
            
    return SA_res   
