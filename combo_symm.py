import argparse
import os
import numpy as np
from prody import *
import math
from math import cos, sin, sqrt
import argparse

def VectorAlign(a, b):
    v = np.cross(a, b)
    c = np.dot(a, b)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    return I + vx + np.dot(vx, vx)/(1 + c)

def RotateMatrix(p1, p2, theta):
    from math import cos, sin, sqrt
    N = (p2 - p1)
    Nm = np.linalg.norm(N)
    n = N / Nm
    c = cos(theta)
    t = (1 - cos(theta))
    s = sin(theta)
    n_lis = n.tolist()
    X = n_lis[0]
    Y = n_lis[1]
    Z = n_lis[2]
    d11 = t * X ** 2 + c
    d12 = t * X * Y - s * Z
    d13 = t * X * Z + s * Y
    d21 = t * X * Y + s * Z
    d22 = t * Y ** 2 + c
    d23 = t * Y * Z - s * X
    d31 = t * X * Z - s * Y
    d32 = t * Y * Z + s * X
    d33 = t * Z ** 2 + c
    M = np.array([[d11, d12, d13],
                  [d21, d22, d23],
                  [d31, d32, d33]])
    return M

def rotate_along_XYZaxis(p0,axis,theta):
    Rx = np.array([[1, 0, 0],
                [0, cos(theta), -sin(theta)],
                [0, sin(theta), cos(theta)]])
    Ry = np.array([[cos(theta), 0, sin(theta)],
            [0, 1, 0],
            [-sin(theta), 0, cos(theta)]])
    Rz = np.array([[cos(theta), -sin(theta), 0],
        [sin(theta), cos(theta), 0],
        [0, 0, 1]])
    flip_z = np.array([[0, 1, 0],
        [1, 0, 0],
        [0, 0, -1]])
    Rm = {'X':Rx, 'Y':Ry, 'Z':Rz,'flip_z':flip_z}
    p = np.dot(Rm[axis], p0.T).T
    return p

def rotate2target(aligned_model0,aligned_model1,combo_model,flip,theta):
    model_0 = parsePDB(aligned_model0)
    model_1 = parsePDB(aligned_model1)
    if flip[0] == '1':
        flipped_model_0_coord = rotate_along_XYZaxis(model_0.getCoords(),axis='flip_z',theta=math.pi)
        model_0.setCoords(flipped_model_0_coord)
    if flip[1] == '1':
        flipped_model_1_coord = rotate_along_XYZaxis(model_1.getCoords(),axis='flip_z',theta=math.pi)
        model_1.setCoords(flipped_model_1_coord)

    Ry = np.array([[cos(theta), 0, sin(theta)],
            [0, 1, 0],
            [-sin(theta), 0, cos(theta)]])

    model_0_coord = model_0.getCoords()
    new_model0_coord = np.dot(Ry, model_0_coord.T).T
    model_0.setCoords(new_model0_coord)
    writePDB(combo_model+'_1.pdb',model_0)
    writePDB(combo_model+'_2.pdb',model_1)

if __name__ == '__main__':
    rules = {'C2+C2':
                {'C2':0,
                },
            }
    parser = argparse.ArgumentParser()
    parser.add_argument('-m1', '--motif1', type = str, \
            help='input motif1')
    parser.add_argument('-symm1', '--symm1', type = str, \
            help='symmetry symbol of input motif1')
    parser.add_argument('-m2', '--motif2', type = str, \
            help='input motif2')
    parser.add_argument('-symm2', '--symm2', type = str, \
            help='symmetry symbol of input motif2')
    parser.add_argument('-s', '--target_symm', type = str, \
            help='target axis to be aligned to')
    parser.add_argument('-f', '--flip_state', type = str, \
            help='flip states of the two motifs, e.g., 01')
    parser.add_argument('-o', '--output_dir', type = str, \
            help='output dir to save the start models')
    
    args = parser.parse_args()
    motif1 = args.motif1
    motif2 = args.motif2
    symm1 = args.symm1
    symm2 = args.symm2
    target_symm = args.target_symm
    flip_state = args.flip_state
    id1 = os.path.split(motif1)[1][:4]
    id2 = os.path.split(motif2)[1][:4]
    output_dir = args.output_dir+'/'+id1+'_'+id2+'_'+target_symm+'/'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    combo_model = os.path.join(output_dir,'_'.join([id1,id2,target_symm,flip_state]))
    if symm1+'+'+symm2 in list(rules.keys()):
        solutions = rules[symm1+'+'+symm2]
        theta = np.deg2rad(solutions[target_symm])
        rotate2target(motif1,motif2,combo_model,flip_state,theta)
    else:
        raise ValueError('No supported target symmetry for this SymmCombo: ', symm1+'+'+symm2)