import os
import shutil
import argparse
import numpy as np
from pymol import cmd
from topoly import homfly, gln, alexander
from topoly.params import Closure

def extract_combined_CA_coord(pdb):
    cmd.delete('all')
    xyz_CA = []
    cmd.load(pdb,'0')
    ch_IDs = cmd.get_chains('0 and not chain X and not chain Y')
    for ch in ch_IDs:
        cmd.select('ch_CA','chain '+ch+' and name CA')
        ch_ca_xyz = cmd.get_coords('ch_CA', '1')

        xyz_CA.append(ch_ca_xyz)
    cmd.delete('all')
    return xyz_CA

def write_xyz(coord,output_xyz):
    with open(output_xyz,'w') as f:
        ch_num = len(coord)
        for i,ch in enumerate(coord):
            last_line = ''
            for atom in ch:
                curr_line = ' '.join([str(k) for k in atom])+'\n'
                if curr_line == last_line:
                    continue
                else:
                    f.write(curr_line)
                    last_line = curr_line
            if i < ch_num-1:
                f.write('X\n')

def topoly_multimer_homfly_closed(topo_type,xyz):
    if topo_type == 'knot':
        closed_chain = homfly(xyz, closure=Closure.CLOSED,chiral = True,max_cross=500)
    elif topo_type == 'link':
        closed_chain = homfly(xyz, closure=Closure.CLOSED,chiral = True,max_cross=500)
    return closed_chain

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sc', '--symm_combo_flip', type = str, \
            help='input symm combo flip')
    args = parser.parse_args()
    symm_combo_flip = args.symm_combo_flip

    target_dir =  './outputs/04_topo/'+symm_combo_flip+'/'
    pdb_list = [k for k in os.listdir(target_dir) if k.startswith(symm_combo_flip) and k.endswith('combined_whole.pdb')]
    result_dir = './outputs/05_result/'+symm_combo_flip+'/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    code_map = {}
    for pdb in pdb_list:
        xyz_CA = extract_combined_CA_coord(target_dir+pdb)
        if len(xyz_CA) == 1:
            topo_type = 'knot'
        elif len(xyz_CA) > 1:
            topo_type = 'link'
        name = pdb[:-4]
        xyz_f = result_dir+name+'_ca.xyz'
        write_xyz(xyz_CA,xyz_f)
        topo_result = topoly_multimer_homfly_closed(topo_type,xyz_f)
        print (name,topo_result)
        if len(topo_result) > 10:
            if topo_result not in list(code_map.keys()):
                curr_num = len(list(code_map.keys()))
                code_map[topo_result] = str(curr_num)+'_'+topo_type
            topo_result = code_map[topo_result]
        else:
            code_map[topo_result] = topo_result
        shutil.copy(target_dir+pdb,result_dir+name+'_'+topo_result+'.pdb')
    with open(result_dir+'code_map.txt','w') as f:
        for code in list(code_map.keys()):
            f.write('{}\t{}\n'.format(code_map[code],code))