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

from Jwalk import PDBTools, GridTools, SATools
import math, os, sys, argparse
import time

global lys_CA, acid_CA, cys_CA, xl_CA, xl_pair, surface
global lysines, acidic, cysteines
global lys, acid, cys, xl_starters, lys_key
global dens_map
global vox_size
global max_dist
    
def find_specific_xl(k_list):
    
    k1 = k_list[0],k_list[1],k_list[2]
    
    some_xl = {}
        
    comb = [[1,0,0],[-1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1],
            [1,0,1],[-1,0,1],[0,1,1],[0,-1,1],[1,-1,0],[-1,-1,0],[1,1,0],[-1,1,0],[1,0,-1],[-1,0,-1],[0,1,-1],[0,-1,-1],
            [1,1,1],[1,-1,1],[-1,1,1],[-1,-1,1],[1,1,-1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]
            ]
            
    far = (math.sqrt((vox_size**2)*2))  
    furthest = (math.sqrt((vox_size**2)*3))
    
    queue = []
    visited = {} # list works as all the coordinates that have been visited - dictionary gives the path to said coordinate from startpoint
    add_length = {}
    
    
    for j in xl_starters[k1]:
        queue.append([j[0],j[1],j[2]])
        visited[j[0],j[1],j[2]] = [[j[0],j[1],j[2]]]
        add_length[j[0],j[1],j[2]] = 0
                    
    while queue:         
        x_n,y_n,z_n = queue.pop(0)
        if add_length[x_n,y_n,z_n] <= max_dist:           
            for c in comb:
                x_temp = x_n+c[0]
                y_temp = y_n+c[1]    
                z_temp = z_n+c[2]
                if (x_temp,y_temp,z_temp) not in visited:           
                    if ((0 <= x_temp < dens_map.x_size()) and (0 <= y_temp < dens_map.y_size()) and (0 <= z_temp < dens_map.z_size())):                                                   
                        temp_list=visited[x_n,y_n,z_n][:]
                        temp_list.append([x_temp,y_temp,z_temp])
                        visited[x_temp,y_temp,z_temp] = temp_list # updated visited list
                        
                        if dens_map.fullMap[z_temp][y_temp][x_temp] <=0:   # if the voxel is in empty space
                            queue.append(([x_temp,y_temp,z_temp]))
                        # calculate the distance
                        diff_x = x_temp - x_n
                        diff_y = y_temp - y_n
                        diff_z = z_temp - z_n
                        if diff_x != 0 and diff_y != 0 and diff_z != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ furthest
                        elif diff_x != 0 and diff_y != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ far
                        elif diff_x != 0 and diff_z != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ far
                        elif diff_y != 0 and diff_z != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ far   
                        else:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ vox_size
    
                        
    # now we have a full set of paths into empty space starting from k1 - now need to extract paths to other lysines.
    
    for aa1 in xl_pair:
        
        h = aa1[0],aa1[1],aa1[2]
        if k1 == h:
            
            k2 = xl_pair[aa1][0],xl_pair[aa1][1],xl_pair[aa1][2]
            best_xl = 9999
            xl_length_dict = {}
            
            for j in xl_starters[k2]:
                
                (x,y,z) = j
                
                if (x,y,z) in visited:
                    
                    visited[(x,y,z)].insert(0,xl_CA[k1])
                    visited[(x,y,z)].append(xl_CA[k2])
                    
                    # now get the lengths   - all will start at k1 and end at k2.
                    
                    xl_length = 0
                        
                    for i in [1,len(visited[(x,y,z)])-1]:
                        (x_1,y_1,z_1) = visited[(x,y,z)][i-1]
                        (x_2,y_2,z_2) = visited[(x,y,z)][i]
                        add_length[(x,y,z)] += math.sqrt((x_1-x_2)**2 + (y_1-y_2)**2 + (z_1-z_2)**2)
                        
                    xl_length = add_length[(x,y,z)]

                    xl_length_dict[xl_length] = visited[(x,y,z)] # dictionary of length:path
                    
                    
                    
                    if best_xl > xl_length:
                        best_xl = xl_length
            
            # now adding shortest xl to the final list
             
            if best_xl != 9999:
                    # this is just to order the dict so that chain goes alphabetically
                some_xl[k1,k2,best_xl] = xl_length_dict[best_xl]     # start lys, end lys, length of xl = path of xl  
    
    return some_xl 
           
def find_xl_2type(k1):
    
    
    some_xl = {} 
       
    comb = [[1,0,0],[-1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1],
            [1,0,1],[-1,0,1],[0,1,1],[0,-1,1],[1,-1,0],[-1,-1,0],[1,1,0],[-1,1,0],[1,0,-1],[-1,0,-1],[0,1,-1],[0,-1,-1],
            [1,1,1],[1,-1,1],[-1,1,1],[-1,-1,1],[1,1,-1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]
            ]
            
    far = (math.sqrt((vox_size**2)*2))  
    furthest = (math.sqrt((vox_size**2)*3))
    
    queue = []
    visited = {} # list works as all the coordinates that have been visited - dictionary gives the path to said coordinate from startpoint
    add_length = {}
    for j in lys[k1]:
        
        queue.append([j[0],j[1],j[2]])
        visited[j[0],j[1],j[2]] = [[j[0],j[1],j[2]]]
        add_length[j[0],j[1],j[2]] = 0
                   
    while queue:         
        x_n,y_n,z_n = queue.pop(0)
        if add_length[x_n,y_n,z_n] <= max_dist:           
            for c in comb:
                x_temp = x_n+c[0]
                y_temp = y_n+c[1]    
                z_temp = z_n+c[2]
                if (x_temp,y_temp,z_temp) not in visited:           
                    if ((0 <= x_temp < dens_map.x_size()) and (0 <= y_temp < dens_map.y_size()) and (0 <= z_temp < dens_map.z_size())):                                                   
                        temp_list=visited[x_n,y_n,z_n][:]
                        temp_list.append([x_temp,y_temp,z_temp])
                        visited[x_temp,y_temp,z_temp] = temp_list # updated visited list
                        
                        if dens_map.fullMap[z_temp][y_temp][x_temp] <=0:   # if the voxel is in empty space
                            queue.append(([x_temp,y_temp,z_temp]))
                        # calculate the distance
                        diff_x = x_temp - x_n
                        diff_y = y_temp - y_n
                        diff_z = z_temp - z_n
                        if diff_x != 0 and diff_y != 0 and diff_z != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ furthest
                        elif diff_x != 0 and diff_y != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ far
                        elif diff_x != 0 and diff_z != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ far
                        elif diff_y != 0 and diff_z != 0:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ far   
                        else:
                            add_length[x_temp,y_temp,z_temp] = add_length[x_n,y_n,z_n]+ vox_size

                        
    # now we have a full set of paths into empty space starting from k1 - now need to extract paths to other lysines.
    for k2 in acid:
        if k1 != k2:
            best_xl = 9999
            xl_length_dict = {}
                
            for j in acid[k2]:    # j = cycling through possible end coords of k2
                
                (x,y,z) = j
                
                if (x,y,z) in visited:
                    
                    visited[(x,y,z)].insert(0,lys_CA[k1])
                    visited[(x,y,z)].append(acid_CA[k2])
                    
                    # now get the lengths   - all will start at k1 and end at k2.
                    
                    xl_length = 0
                        
                    for i in [1,len(visited[(x,y,z)])-1]:
                        (x_1,y_1,z_1) = visited[(x,y,z)][i-1]
                        (x_2,y_2,z_2) = visited[(x,y,z)][i]
                        add_length[(x,y,z)] += math.sqrt((x_1-x_2)**2 + (y_1-y_2)**2 + (z_1-z_2)**2)
                        
                    xl_length = add_length[(x,y,z)]
    
                    xl_length_dict[xl_length] = visited[(x,y,z)] # dictionary of length:path
                    
                    if best_xl > xl_length:
                        best_xl = xl_length
            
            # now adding shortest xl to the final list
                
            if best_xl != 9999:
                    if k1[1] < k2[1]:             # this is just to order the dict so that chain goes alphabetically
                        some_xl[k1,k2,best_xl] = xl_length_dict[best_xl]     # start lys, end lys, length of xl = path of xl  
                    elif k2[1] < k1[1]:
                        some_xl[k2,k1,best_xl] = xl_length_dict[best_xl]
                    elif k1[0] <  k2[0]:  
                        some_xl[k1,k2,best_xl] = xl_length_dict[best_xl]
                    else:
                        some_xl[k2,k1,best_xl] = xl_length_dict[best_xl] 
            
    return some_xl      
                              
def paralell_BFS(uv_list = False):
    import multiprocessing
    from multiprocessing import Pool,  Process, freeze_support
    import itertools
        
    freeze_support()
    
    if uv_list == True:
        final_XL = {}       
        pool = Pool(processes=multiprocessing.cpu_count()) 
        xl_dictionaries = pool.map(find_specific_xl, xl_pair)
        for c in xl_dictionaries:
            final_XL.update(c)
        
    
    else:
        final_XL = {}
        pool = Pool(processes=multiprocessing.cpu_count()) 
        xl_dictionaries = pool.map(find_xl_2type, lys)
        for c in xl_dictionaries:
            final_XL.update(c)        
    
    return final_XL

################################

lysines = False 
types2 = False 
uv_list = False
max_dist = 60
vox = 1
surface = False

parser = argparse.ArgumentParser(description='JWALK: Calculate SASDs on your target PDB files')

parser.add_argument('-lys', action= "store_true",
               help='calculate lysine crosslinks (default)')
parser.add_argument('-xl_list', nargs = 1,
               help='calculate crosslinks from input list')
parser.add_argument('-i',  nargs=1,
               help='specify input pdb: -i <inputfile.pdb>')
parser.add_argument('-aa1',  nargs=1,
               help='specify start amino acid (three letter code e.g. LYS)')
parser.add_argument('-aa2',  nargs=1,
               help='specify end amino acid (three letter code e.g. LYS)')
parser.add_argument('-max_dist',  nargs=1,
               help='specify maximum crosslink distance in Angstroms')
parser.add_argument('-vox',  nargs=1,
               help='specify voxel size of grid')
parser.add_argument('-surface',  action = "store_true",
               help='use higher accuracy method to calculate solvent accessibility - slower')

args = parser.parse_args()
if args.lys:
    lysines = True
    aa_1 = "LYS" ; aa_2 = "LYS"
    types2 = True
    
elif args.xl_list:
    uv_list = True
    uv_xl = args.xl_list[0]

elif args.aa1:
    if args.aa2:
        types2 = True
        aa_1 = args.aa1[0].upper()
        aa_2 = args.aa2[0].upper()
        
    else:
        print "Please specify both aa1 AND aa2 if you want to use this option"
        sys.exit(2)
    
else:
    lysines = True
    aa_1 = "LYS" ; aa_2 = "LYS"
    types2 = True
    
if args.max_dist:
    max_dist = int(args.max_dist[0])
        
if args.i:
    pdb_list = [args.i[0]]
else:
    pdb_list = [i for i in os.listdir('./') if i.endswith('.pdb')]

if args.vox:
    vox = int(args.vox[0])
    
if args.surface:
    surface = True

#######################################

for pdb in pdb_list:

    print "calculating crosslinks on", pdb
    
    structure_instance=PDBTools.read_PDB_file(pdb)
    
    grid = GridTools.makeGrid(structure_instance,vox,3)
    
    if uv_list == True:
        
        xl_pair, xl_CA = GridTools.mark_CA_UV(grid, structure_instance,uv_xl)
        
        if surface == True:
            
            xl_CA = SATools.calculate_if_SA(structure_instance,xl_CA,uv_list)
            dens_map = GridTools.generate_SAS_UV(grid,structure_instance,xl_CA)
            xl_starters = GridTools.find_surface_voxels2(xl_CA,dens_map)
            
        else:
            dens_map = GridTools.generate_SAS_UV(grid,structure_instance,xl_CA)
            
            xl_starters,rem_x = GridTools.find_surface_voxels_surface_false(xl_CA,dens_map,uv_list)
            
            xl_pair = GridTools.remove_crosslinks(xl_pair,rem_x)
            
            
        vox_size = dens_map.apix

        crosslinks = paralell_BFS(uv_list)
        
        crosslinks_final = GridTools.remove_duplicates(crosslinks,uv_list)
        
        PDBTools.write_xl_to_txt(crosslinks_final,pdb)
        PDBTools.write_xl_to_pdb(dens_map,crosslinks_final,pdb)
        print len(crosslinks_final), "SASDs calculated"
    
    elif types2 == True:
        lys_CA, acid_CA = GridTools.mark_CA_2types(grid, structure_instance,aa_1,aa_2)
        if surface == True:
            lys_CA  = SATools.calculate_if_SA(structure_instance,lys_CA)
            acid_CA = SATools.calculate_if_SA(structure_instance,acid_CA)
            dens_map = GridTools.generate_SAS_2type(grid,structure_instance,lys_CA,acid_CA)
            lys,acid = GridTools.find_surface_voxels_2type(lys_CA,acid_CA,dens_map)
        else:
            dens_map = GridTools.generate_SAS_2type(grid,structure_instance,lys_CA,acid_CA)
            lys,rem_x = GridTools.find_surface_voxels_surface_false(lys_CA,dens_map)
            acid,rem_x = GridTools.find_surface_voxels_surface_false(acid_CA,dens_map)
        
           
        vox_size = dens_map.apix
        
        crosslinks = paralell_BFS()
        
        crosslinks_final = GridTools.remove_duplicates(crosslinks)
        
        PDBTools.write_xl_to_txt(crosslinks_final,pdb)
        PDBTools.write_xl_to_pdb(dens_map,crosslinks_final,pdb)
        print len(crosslinks_final), "SASDs calculated"
    
    


