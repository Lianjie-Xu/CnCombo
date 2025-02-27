import argparse
import os
import subprocess
import time

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input_pdb1', type = str, \
            help='input multimer pdb file 1')
    parser.add_argument('-i2', '--input_pdb2', type = str, \
            help='input multimer pdb file 2')
    parser.add_argument('-t', '--target_axis', type = str, default='0,0,1', \
            help='target axis to be aligned to')
    
    parser.add_argument('-symm1', '--symm1', type = str, \
            help='symmetry symbol of input motif1')
    parser.add_argument('-symm2', '--symm2', type = str, \
            help='symmetry symbol of input motif2')
    parser.add_argument('-s', '--target_symm', type = str, \
            help='target axis to be aligned to')
    parser.add_argument('-c', '--cutoff', type = float, default=5.0, \
            help='cutoff for translation')
    parser.add_argument('-g', '--interval', type = float, default=5.0, \
            help='interval for rotation')
    args = parser.parse_args()

    cutoff = args.cutoff
    interval = args.interval
    input_pdb1 = args.input_pdb1
    input_pdb2 = args.input_pdb2
    axis_ch1,axis_ch2 = 'X','Y'
    motif_dir = './motifs/'
    id1,id2 = os.path.split(input_pdb1)[-1][:4],os.path.split(input_pdb2)[-1][:4]
    motif1 = motif_dir+id1+'_Zaxis.pdb'
    motif2 = motif_dir+id2+'_Zaxis.pdb'
    p1_1 = subprocess.Popen(f'python combo_multimer.py --input_pdb {input_pdb1} --axis_ch {axis_ch1} --output_pdb {motif1}', shell=True)
    p1_1.wait()
    p1_2 = subprocess.Popen(f'python combo_multimer.py --input_pdb {input_pdb2} --axis_ch {axis_ch2} --output_pdb {motif2}', shell=True)
    p1_2.wait()

    symm1 = args.symm1
    symm2 = args.symm2
    target_symm = args.target_symm
    assert symm1 == symm2 == target_symm and symm1[0] == symm2[0] == target_symm[0] == 'C', "Currently only Cn+Cn=Cn is supported..."
    n = target_symm[1:]
    flip_states = ['00','01','10','11']
    start_dir = './outputs/01_start/'
    if not os.path.exists(start_dir):
            os.makedirs(start_dir)
    for flip in flip_states:
        print ('+'*50+'Flip state: {}'.format(flip)+'+'*50)
        p2 = subprocess.Popen('python combo_symm.py -m1 {0} -symm1 {1} -m2 {2} -symm2 {3} -s {4} -f {5} -o {6}'.format(motif1,symm1,motif2,symm2,target_symm,flip,start_dir), shell=True)
        p2.wait()
        symm_combo_flip = '_'.join([id1,id2,target_symm,flip])
        symm_combo = '_'.join([id1,id2,target_symm])
        start_sub_dir = start_dir+symm_combo+'/'
        start1 = start_sub_dir+'_'.join([id1,id2,target_symm,flip,'1'])+'.pdb'
        start2 = start_sub_dir+'_'.join([id1,id2,target_symm,flip,'2'])+'.pdb'
        trans_rot_dir = './outputs/02_trans_rot/'
        if not os.path.exists(trans_rot_dir):
            os.mkdir(trans_rot_dir)
        p3 = subprocess.Popen('python combo_trans_rot_Cn.py -n {0} -m1 {1} -m2 {2} -g {3} -c {4} -o {5}'.format(n,start1,start2,interval,cutoff,trans_rot_dir), shell=True)
        p3.wait()

        os.chdir(trans_rot_dir+symm_combo_flip)
        p4 = subprocess.Popen('python ../../../combo_linkages_Cn.py -sc {}'.format(symm_combo_flip), shell=True)
        p4.wait()

        os.chdir(os.path.abspath(os.path.join(os.getcwd(), "../../..")))
        p5 = subprocess.Popen('python combo_topo_Cn.py -sc {}'.format(symm_combo_flip), shell=True)
        p5.wait()
    
    end = time.time()
    print ('Running time: %.2f minutes' % ((end-start)/60))