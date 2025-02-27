from itertools import tee
import os
import numpy as np
from prody import *
from pymol import cmd
import math
from math import cos, sin, sqrt
import argparse
from pyrosetta import *
init('-mute all')

def get_term_index(pdb):
    term_index_dict = {}
    pose = pose_from_pdb(pdb)
    chains = list(pyrosetta.rosetta.core.pose.get_chains(pose))
    for i,chain in enumerate(chains):
        sub_pose = pyrosetta.rosetta.core.pose.Pose.split_by_chain(pose,chain)
        ch_id = sub_pose.pdb_info().chain(1)
        total = sub_pose.total_residue()
        start = sub_pose.pdb_info().number(1)
        end = sub_pose.pdb_info().number(total)
        term_index_dict[ch_id] = [start,end]
    return term_index_dict

def angle_between(p1, p2):
    ang1 = np.arctan2(*p1[::-1])
    ang2 = np.arctan2(*p2[::-1])
    return np.rad2deg((ang1 - ang2) % (2 * np.pi))

def distance_between(p1, p2):
    return np.linalg.norm(p1-p2)

def translate2cutoff_Cn(n,aligned_model0,aligned_model1,fusion_name,fused_model_dir,interval,cutoff=5.0):
    cmd.delete('all')
    cmd.load(aligned_model0,'0')
    cmd.load(aligned_model1,'1')
    model0_chains = cmd.get_chains('0 and not chain X and not chain Y')
    for i,ch in enumerate(model0_chains):
        cmd.select('0_subunit_{}'.format(i),'0 and chain {}'.format(ch))
    for i,ch in enumerate(model0_chains):
        cmd.alter('0_subunit_{}'.format(i),'chain="{}"'.format(chr(i+97)))
    model1_chains = cmd.get_chains('1 and not chain X and not chain Y')
    for i,ch in enumerate(model1_chains):
        cmd.select('1_subunit_{}'.format(i),'1 and chain {}'.format(ch))
    for i,ch in enumerate(model1_chains):
        cmd.alter('1_subunit_{}'.format(i),'chain="{}"'.format(chr(i+65)))
    
    cmd.select('model_0_ca','0 and not chain X and not chain Y and name CA')
    smallest_z_model0 = np.min(cmd.get_coords('model_0_ca', 1), axis=0)[2]
    cmd.select('model_1_ca','1 and not chain X and not chain Y and name CA')
    biggest_z_model1 = np.max(cmd.get_coords('model_1_ca', 1), axis=0)[2]
    trans_d = biggest_z_model1 - smallest_z_model0
    cmd.translate([0,0,trans_d+cutoff],'0')
    
    model0_term_index_dict = get_term_index(aligned_model0)
    model0_NT, model0_CT = model0_term_index_dict['A']
    model1_term_index_dict = get_term_index(aligned_model1)
    model1_NT, model1_CT = model1_term_index_dict['A']
    
    cmd.select('0-CT','0 and chain {} and resi {} and name CA'.format('a', model0_CT))
    cmd.select('1-NT','1 and chain {} and resi {} and name CA'.format('A', model1_NT))
    m0_CT_coord = cmd.get_coords('0-CT', 1)[0]
    m1_NT_coord = cmd.get_coords('1-NT', 1)[0]
    m0_CT_coord_xy = m0_CT_coord[:2]
    m1_NT_coord_xy = m1_NT_coord[:2]
    rot_ang0_m0_m1 = angle_between(m1_NT_coord_xy, m0_CT_coord_xy)
    cmd.copy_to('1_start_m0_m1','1')
    cmd.rotate('z', -rot_ang0_m0_m1, '1_start_m0_m1',origin=[0,0,0])
    m1_NT_new_coord = cmd.get_coords('1_start_m0_m1 and chain {} and resi {} and name CA'.format('A',model1_NT), 1)[0]
    minD_m0_CT_m1_NT = distance_between(m0_CT_coord, m1_NT_new_coord)
    
    cmd.select('1-CT','1 and chain {} and resi {} and name CA'.format('A',model1_CT))
    cmd.select('0-NT','0 and chain {} and resi {} and name CA'.format('a',model0_NT))
    m1_CT_coord = cmd.get_coords('1-CT', 1)[0]
    m0_NT_coord = cmd.get_coords('0-NT', 1)[0]
    m1_CT_coord_xy = m1_CT_coord[:2]
    m0_NT_coord_xy = m0_NT_coord[:2]
    rot_ang0_m1_m0 = angle_between(m1_CT_coord_xy, m0_NT_coord_xy)
    cmd.copy_to('1_start_m1_m0','1')
    cmd.rotate('z', -rot_ang0_m1_m0, '1_start_m1_m0',origin=[0,0,0])
    m1_CT_new_coord = cmd.get_coords('1_start_m1_m0 and chain {} and resi {} and name CA'.format('A',model1_CT), 1)[0]
    minD_m1_CT_m0_NT = distance_between(m0_NT_coord, m1_CT_new_coord)

    if minD_m0_CT_m1_NT < minD_m1_CT_m0_NT:
        m1_start_final = '1_start_m0_m1'
    elif minD_m0_CT_m1_NT > minD_m1_CT_m0_NT:
        m1_start_final = '1_start_m1_m0'
    else:
        raise ValueError('Equal distance between m0->m1 and m1->m0')
    print (m1_start_final, minD_m0_CT_m1_NT, minD_m1_CT_m0_NT)
    cmd.save(os.path.join(os.path.join(fused_model_dir,'../'), '{}_start.pdb'.format(fusion_name)), '0 or {}'.format(m1_start_final))
    
    angle_val = np.linspace(-180/n, 180/n, int(360/(n*interval)+1), endpoint=True)
    print (angle_val)
    id_number = 1
    for i in np.arange(len(angle_val)):
        cmd.copy_to('1_copy_{}'.format(i),m1_start_final)
        cmd.rotate('z', angle_val[i], '1_copy_{}'.format(i),origin=[0,0,0])
        cmd.save(os.path.join(fused_model_dir, '{}_angle_{}_{}.pdb'.format(fusion_name,int(angle_val[i]),id_number)), '0 or 1_copy_{}'.format(i))
        id_number += 1
    cmd.delete('all')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--Cn', type = str, \
            help='Cn+Cn=Cn')
    parser.add_argument('-m1', '--motif1', type = str, \
            help='input motif1')
    parser.add_argument('-m2', '--motif2', type = str, \
            help='input motif2')
    parser.add_argument('-c', '--cutoff', type = float, default=5.0, \
            help='cutoff for translation')
    parser.add_argument('-g', '--interval', type = float, default=5.0, \
            help='interval for rotation')
    parser.add_argument('-o', '--output_dir', type = str, \
            help='output dir to save the fusion model')
    args = parser.parse_args()
    n = int(args.Cn)
    aligned_model0 = args.motif1
    aligned_model1 = args.motif2
    cutoff = args.cutoff
    interval = args.interval
    fusion_name = '_'.join(os.path.split(aligned_model0)[-1].split('_')[:4])
    fused_model_dir = os.path.join(args.output_dir,fusion_name)
    if not os.path.exists(fused_model_dir):
        os.mkdir(fused_model_dir)
    translate2cutoff_Cn(n,aligned_model0,aligned_model1,fusion_name,fused_model_dir,interval,cutoff)
