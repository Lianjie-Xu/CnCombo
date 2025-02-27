import os
import re
import shutil
import argparse
import numpy as np
from prody import *
import math
from pymol import cmd
from pymol import stored
from chempy import cpv
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
init('-mute all -detect_disulf 0')

def clean_pdbx(pdb):
    pose = pose_from_pdb(pdb)
    pose.dump_pdb(os.path.split(pdb)[0]+os.path.split(pdb)[1][:-4]+'_dump.pdb')


def centroid(selection='all', center=0, quiet=1):
    # https://wiki.pymol.org/index.php/Centroid
    model = cmd.get_model(selection)
    nAtom = len(model.atom)

    centroid = cpv.get_null()

    for a in model.atom:
        centroid = cpv.add(centroid, a.coord)
    centroid = cpv.scale(centroid, 1. / nAtom)

    if not int(quiet):
        print(' centroid: [%8.3f,%8.3f,%8.3f]' % tuple(centroid))

    if int(center):
        cmd.alter_state(1, selection, "(x,y,z)=sub((x,y,z), centroid)",
                        space={'centroid': centroid, 'sub': cpv.sub})

    return centroid


def get_axis(multimer):
    result = os.popen("./helper_scripts/symd1.61-linux " + multimer).read()
    return result

def add_axis(old_pdb,trfm_pdb,axis_ch,added_pdb):
    with open(added_pdb,'w') as f:
        with open(old_pdb, 'r') as old:
            old_lines = old.readlines()
            for line in old_lines:
                if (line.startswith('ATOM') or line.startswith('HETATM')) and len(line)>=73:
                    new_line = line[:72]+' '*4+line[76:]
                    f.write(new_line)
                elif line.startswith('END'):
                    continue
                else:
                    f.write(line)
            print ('*'*20)
            with open(trfm_pdb,'r') as trfm:
                trfm_lines = trfm.readlines()
                axis_xyz = []
                for i in range(-8, -5, 1):
                    print (trfm_lines[i])
                    f.write(trfm_lines[i][:21]+axis_ch+trfm_lines[i][22:])
                    axis_xyz.append([float(k) for k in trfm_lines[i].split()[-3:]])
                for i in range(-5, 0, 1):
                    print (trfm_lines[i])
                    f.write(trfm_lines[i])
    print ('Axis: ',axis_xyz)
    return axis_xyz

def VectorAlign(a, b):
    v = np.cross(a, b)
    c = np.dot(a, b)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    return I + vx + np.dot(vx, vx)/(1 + c)

def is_parallel(vec1, vec2):
    assert isinstance(vec1, np.ndarray)
    assert isinstance(vec2, np.ndarray)
    assert vec1.shape == vec2.shape

    vec1_normalized = vec1 / np.linalg.norm(vec1)
    vec2_normalized = vec2 / np.linalg.norm(vec2)

    if 1.0 - abs(np.dot(vec1_normalized, vec2_normalized)) < 1e-6:
        return True
    else:
        return False

def align_symm_axis(pdb,axis_chain,target_axis,output_pdb):
    multimer = parsePDB(pdb)
    axis_coords = multimer.select(' '.join(['chain', axis_chain, 'ca', 'resnum', '1 : 4'])).getCoords()
    vector = axis_coords[2] - axis_coords[0]

    rotate2z = VectorAlign(vector/np.linalg.norm(vector),target_axis)
    struc_Coords = np.dot(rotate2z,(multimer.getCoords() - axis_coords[0]).T).T
    multimer.setCoords(struc_Coords)
    writePDB(output_pdb,multimer)

def CM2origin(pdb,output_pdb):
    cmd.delete('all')
    cmd.load(pdb,'0')
    centroid('0',center=1)
    print ('centroid to origin...')
    cmd.save(output_pdb)
    cmd.delete('all')    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_pdb', type = str, \
            help='input multimer pdb file')
    parser.add_argument('-a', '--axis_ch', type = str, \
            help='chain ID of axis chain to be added into pdb file')
    parser.add_argument('-t', '--target_axis', type = str, default='0,0,1', \
            help='target axis to be aligned to')
    parser.add_argument('-o', '--output_pdb', type = str, \
            help='output pdb to be saved.')
    
    args = parser.parse_args()
    input_pdb = args.input_pdb
    axis_ch = args.axis_ch
    target_axis =  args.target_axis
    output_pdb = args.output_pdb

    get_axis(input_pdb)

    pdb_old = input_pdb
    pdb_trfm = os.path.split(input_pdb)[-1][:-4] + '-trfm.pdb'
    pdb_added = os.path.split(input_pdb)[-1][:-4] + '_AxisAdded.pdb'
    axis_xyz = add_axis(pdb_old,pdb_trfm,axis_ch,pdb_added)
    axis_vec = np.array(axis_xyz[0])-np.array(axis_xyz[-1])

    t_axis_xyz = target_axis.split(',')
    if len(t_axis_xyz) == 3:
        t_axis = [float(k) for k in target_axis.split(',')]
        if not is_parallel(axis_vec, np.array(t_axis)):
            align_symm_axis(pdb_added,axis_ch,t_axis,output_pdb)
            CM2origin(output_pdb,output_pdb)
        else:
            print ('symm axis already aligned...')
            CM2origin(pdb_added,output_pdb)
    else:
        raise ValueError("Unsopported format of target axis, example input: '0,0,1'...")
    
    motif_dir = os.path.dirname(input_pdb)
    shutil.move(pdb_trfm,motif_dir+'/'+pdb_trfm)
    shutil.move(pdb_added,motif_dir+'/'+pdb_added)
    shutil.move(os.path.split(input_pdb)[-1][:-4] + '-info.txt', motif_dir+'/'+os.path.split(input_pdb)[-1][:-4] + '-info.txt')
    shutil.move(os.path.split(input_pdb)[-1][:-4] + '-best.fasta', motif_dir+'/'+os.path.split(input_pdb)[-1][:-4] + '-best.fasta')