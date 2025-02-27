import os
import shutil
import time
import argparse
from prody import *
from pymol import cmd
from pyrosetta import *
from pyrosetta.toolbox import *
init('-mute all')

def write_xl(n,pdb,combo,xl_dir):
    term_index_dict = {}
    pose = pose_from_pdb(pdb)
    chains = list(pyrosetta.rosetta.core.pose.get_chains(pose))
    for i,chain in enumerate(chains):
        sub_pose = pyrosetta.rosetta.core.pose.Pose.split_by_chain(pose,chain)
        ch_id = sub_pose.pdb_info().chain(1)
        total = sub_pose.total_residue()
        start = sub_pose.pdb_info().number(1)
        end = sub_pose.pdb_info().number(total)
        term_index_dict[ch_id] = [str(start),str(end)]
    cha_NT,cha_CT = term_index_dict['a']
    with open(xl_dir+combo+'_xl.txt','w') as xl_txt:
        for i in range(n):
            ch = chr(65+i)
            ch_NT,ch_CT = term_index_dict[ch]
            a_ch_xl = '|'.join([cha_CT,'a',ch_NT,ch])
            xl_txt.write(a_ch_xl+'\n')
            ch_a_xl = '|'.join([ch_CT,ch,cha_NT,'a'])
            xl_txt.write(ch_a_xl+'\n')

def Jwalk_run(jwalk_path,pdb,xl_txt_dir,xl_txt,target_dir):
    xl_txt_f = xl_txt_dir+xl_txt
    run_output = os.popen('python {} -max_dist 1000 -xl_list {} -i {}'.format(jwalk_path,xl_txt_f,pdb)).read()  
    print (run_output)
    # time.sleep(10)

def parse_crosslink(n,name,xl_result_txt,xl_result_pdb,xl_temp_pdb,combo_pdb,altered_pdb,xl_dir,target_dir):
    xl_result = []
    with open(xl_result_txt,'r') as f:
        lines = f.read().strip().split('\n')[1:]
        for i,line in enumerate(lines):
            toks = line.split()
            ch1 = toks[2].split('-')[2]
            ch2 = toks[3].split('-')[2]
            res1 = toks[2].split('-')[1]
            res2 = toks[3].split('-')[1]
            sasd = float(toks[-1])
            xl_result.append([ch1+'-'+ch2,res1,res2,sasd])
    xl_min = min(xl_result,key=lambda x:x[-1])
    xl_min_index = xl_result.index(xl_min)
    print ('The first set of SASD (fusion linker) are: ', xl_result)

    with open(xl_temp_pdb, 'w') as tmp:
        tmp_ch = '0'
        with open(xl_result_pdb,'r') as xl_pdb:
            xl_pdb_data = xl_pdb.read().strip('END\n').split('END')
            xl_ch = xl_pdb_data[xl_min_index]
            xl_ch_lines = [l for l in xl_ch.strip().split('\n') if l[:4]=='ATOM']
            for line in xl_ch_lines[1:-1]:
                new_line = line[:14]+'A'+line[15:21]+tmp_ch+line[22:]
                tmp.write(new_line+'\n')
            xl_number = str(int(xl_ch_lines[-1][22:26])-2)
            tmp.write('TER\n')

    cmd.delete('all')
    cmd.load(combo_pdb,'0')
    cmd.load(xl_temp_pdb,'1')
    cmd.alter('1 and chain '+tmp_ch,'resi=str(int(resi)-1)')
    min_linkage,min_res1,min_res2,min_sasd = xl_min
    C_,N_ = min_linkage.split('-')

    cmd.alter('1 and chain '+tmp_ch,'resi=str(int(resi)+1-1'+'+'+min_res1+')')
    cmd.alter('1 and chain '+tmp_ch,'chain="'+C_+'"')
    cmd.alter('0 and chain '+N_,'resi=str(int(resi)+1-'+min_res2+'+'+xl_number+'+'+min_res1+')')
    cmd.alter('0 and chain '+N_,'chain="'+C_+'"')

    symm_anchor = {}
    for i in range(n):
        symm_anchor[chr(65+i)] = [chr(65+k) for k in range(n) if k != i]
        symm_anchor[chr(97+i)] = [chr(97+k) for k in range(n) if k != i]
    cmd.copy_to('xl_'.format(i),'0 and chain {} or 1'.format(C_))
    for i in range(n-1):
        cmd.copy_to('xl_{}'.format(i),'0 and chain {} or 1'.format(C_))
        cmd.rotate('z', (i+1)*2*180/n, 'xl_{}'.format(i),origin=[0,0,0])
        cmd.alter('xl_{}'.format(i),'chain="{}"'.format(symm_anchor[C_][i]))
    
    cmd.select('ch2A','chain '+C_)
    for i in range(n-1):
        cmd.select('ch2_{}'.format(i),'chain '+symm_anchor[C_][i])
    for i in range(n-1):
        cmd.alter('ch2_{}'.format(i),'chain="{}"'.format(chr(65+i+1)))
    cmd.alter('ch2A','chain="A"')
    
    cmd.remove('0 and not chain X and not chain Y')
    cmd.remove('1 and not chain X and not chain Y')

    cmd.save(altered_pdb,'all')
    cmd.delete('all')
    cmd.load(altered_pdb,'2')
    cmd.save(altered_pdb,'all')
    cmd.delete('all')

    term_index_dict = {}
    pose = pose_from_pdb(altered_pdb)
    chains = list(pyrosetta.rosetta.core.pose.get_chains(pose))
    for i,chain in enumerate(chains):
        sub_pose = pyrosetta.rosetta.core.pose.Pose.split_by_chain(pose,chain)
        ch_id = sub_pose.pdb_info().chain(1)
        total = sub_pose.total_residue()
        start = sub_pose.pdb_info().number(1)
        end = sub_pose.pdb_info().number(total)
        term_index_dict[ch_id] = [str(start),str(end)]
    chA_NT,chA_CT = term_index_dict['A']
    with open(xl_dir+name+'_second_xl.txt','w') as xl_txt:
        for i in range(n):
            ch = chr(65+i)
            ch_NT,ch_CT = term_index_dict[ch]
            A_ch_xl = '|'.join([chA_CT,'A',ch_NT,ch])
            xl_txt.write(A_ch_xl+'\n')

    Jwalk_run(jwalk_path,os.path.split(altered_pdb)[-1],xl_dir,xl_dir+name+'_second_xl.txt',target_dir)

def parse_second_crosslink(n,xl_second_result_txt,xl_second_result_pdb,xl_temp_pdb,combo_semi_pdb,altered_semi_pdb):
    xl_result = []
    with open(xl_second_result_txt,'r') as f:
        lines = f.read().strip().split('\n')[1:]
        for i,line in enumerate(lines):
            toks = line.split()
            ch1 = toks[2].split('-')[2]
            ch2 = toks[3].split('-')[2]
            res1 = toks[2].split('-')[1]
            res2 = toks[3].split('-')[1]
            sasd = float(toks[-1])
            xl_result.append([ch1+'-'+ch2,res1,res2,sasd])
    xl_min = min(xl_result,key=lambda x:x[-1])
    xl_min_index = xl_result.index(xl_min)
    print ('The second set of SASD (ligation linker) are: ', xl_result)

    with open(xl_temp_pdb, 'w') as tmp:
        tmp_ch = '0'
        with open(xl_second_result_pdb,'r') as xl_pdb:
            xl_pdb_data = xl_pdb.read().strip('END\n').split('END')
            xl_ch = xl_pdb_data[xl_min_index]
            xl_ch_lines = xl_ch.strip().split('\n')  
            for line in xl_ch_lines:
                if line[:4] == 'ATOM':
                    new_line = line[:14]+'A'+line[15:21]+tmp_ch+line[22:]
                    tmp.write(new_line+'\n')
                else:
                    continue
            xl_number = str(int(xl_ch_lines[-1][22:26]))
            tmp.write('TER\n')

    cmd.delete('all')
    cmd.load(combo_semi_pdb,'0')
    cmd.load(xl_temp_pdb,'1')
    min_linkage,min_res1,min_res2,min_sasd = xl_min
    C_,N_ = min_linkage.split('-')
    min_res1,min_res2 = str(int(min_res1)),str(int(min_res2))
    symm_anchor = {}
    for i in range(n):
        symm_anchor[chr(65+i)] = [chr(65+k) for k in range(n) if k != i]
        symm_anchor[chr(97+i)] = [chr(97+k) for k in range(n) if k != i]
    if C_ == N_:
        print ('This is a LINK!')
        cmd.alter('1 and chain '+tmp_ch,'resi=str(int(resi)+1-1'+'+'+min_res1+')')
        cmd.alter('1 and chain '+tmp_ch,'chain="'+C_+'"')
        cmd.copy_to('xl_','0 and chain {} or 1'.format(C_))
        for i in range(n-1):
            cmd.copy_to('xl_{}'.format(i),'0 and chain {} or 1'.format(C_))
            cmd.rotate('z', (i+1)*2*180/n, 'xl_{}'.format(i),origin=[0,0,0])
            cmd.alter('xl_{}'.format(i),'chain="{}"'.format(symm_anchor[C_][i]))
        
        cmd.select('ch2A','chain '+C_)
        cmd.alter('ch2A','chain="A"')
        for i in range(n-1):
            cmd.select('ch2BZ','chain '+symm_anchor[C_][i])
            cmd.alter('ch2BZ','chain="{}"'.format(chr(65+i+1)))
        
        cmd.remove('0 and not chain X and not chain Y')
        cmd.remove('1 and not chain X and not chain Y')

        cmd.save(altered_semi_pdb,'all')
        cmd.delete('all')
        cmd.load(altered_semi_pdb,'2')
        cmd.save(altered_semi_pdb,'all')
        cmd.delete('all')
    else:
        cmd.copy_to('xl_','1 and chain '+tmp_ch)
        cmd.alter('1 and chain '+tmp_ch,'resi=str(int(resi)+1-1'+'+'+min_res1+')')
        cmd.alter('1 and chain '+tmp_ch,'chain="'+C_+'"')

        mono_len = int(min_res1)-int(min_res2)+1
        curr_total = int(xl_number)+mono_len
        if ord(N_) > ord(C_):
            period = ord(N_)-ord(C_)
        else:
            period = n+ord(N_)-ord(C_)
        print ('period: ', period)
        ABC = [chr(65+c) for c in range(n)]*2
        period_repeat = ABC[ord(C_)-65:ord(C_)-65+n]*n
        print ('period_repeat: ', period_repeat)
        if period == 1 or n % period != 0:
            print ('This is a KNOT!')
            for i in range(n-1):
                next_mono_id = period_repeat[period*(i+1)]
                cmd.alter('0 and chain '+next_mono_id,'resi=str(int(resi)+1-'+min_res2+'+'+str(curr_total)+')')
                cmd.alter('0 and chain '+next_mono_id,'chain="'+C_+'"')
                curr_total += mono_len
                cmd.copy_to('xl_{}'.format(i),'xl_')
                cmd.rotate('z', period*(i+1)*2*180/n, 'xl_{}'.format(i), origin=[0,0,0])
                cmd.alter('xl_{}'.format(i),'resi=str(int(resi)+1-1+'+str(curr_total)+')')
                curr_total += int(xl_number)
                cmd.alter('xl_{}'.format(i),'chain="'+C_+'"')

            cmd.select('chA','chain '+C_)
            cmd.alter('chA','chain="A"')
        elif period != 1 and n % period == 0:
            print ('This is a LINK!')
            link_ch_num = n % period
            for i in range(n/link_ch_num-1):
                next_mono_id = period_repeat[period*(i+1)]
                cmd.alter('0 and chain '+next_mono_id,'resi=str(int(resi)+1-'+min_res2+'+'+str(curr_total)+')')
                cmd.alter('0 and chain '+next_mono_id,'chain="'+C_+'"')
                curr_total += mono_len
                cmd.copy_to('xl_{}'.format(i),'xl_')
                cmd.rotate('z', period*(i+1)*2*180/n, 'xl_{}'.format(i), origin=[0,0,0])
                cmd.alter('xl_{}'.format(i),'resi=str(int(resi)+1-1+'+str(curr_total)+')')
                curr_total += int(xl_number)
                cmd.alter('xl_{}'.format(i),'chain="'+C_+'"')
            cmd.copy_to('link_','0 and chain {} or 1'.format(C_))

            for j in range(link_ch_num):
                cmd.copy_to('link_{}'.format(j),'0 and chain {}'.format(C_))
                cmd.rotate('z', (j+1)*2*180/n, 'link_{}'.format(j),origin=[0,0,0])
                cmd.alter('link_{}'.format(j),'chain="{}"'.format(symm_anchor[C_][j]))
            # renamed chain IDs to ABCD...
            cmd.select('chA','chain '+C_)
            cmd.alter('chA','chain="A"')
            for i in range(n-1):
                cmd.select('ch2BZ','chain '+symm_anchor[C_][i])
                cmd.alter('ch2BZ','chain="{}"'.format(chr(65+i+1)))
            # remove extra chains
            cmd.remove('0 and not chain X and not chain Y')
            cmd.remove('1 and not chain X and not chain Y')
        cmd.remove('xl_')
        # save 
        cmd.save(altered_semi_pdb,'all')
        cmd.delete('all')
        cmd.load(altered_semi_pdb,'2')
        cmd.save(altered_semi_pdb,'all')
        cmd.delete('all')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sc', '--symm_combo_flip', type = str, \
            help='input symm combo flip')
    args = parser.parse_args()
    symm_combo_flip = args.symm_combo_flip
    n = int(symm_combo_flip.split('_')[-2].strip('C'))

    jwalk_path = '../../../helper_scripts/Jwalk/Jwalk.v1.1.py'
    pdb_dir = '../../02_trans_rot/'+symm_combo_flip+'/'
    pdb_list = [k for k in os.listdir(pdb_dir) if k.startswith(symm_combo_flip)]

    xl_dir = '../../03_crosslink/'+symm_combo_flip+'/'
    if not os.path.exists(xl_dir):
        os.makedirs(xl_dir)
    xl_txt = '_'.join(symm_combo_flip.split('_')[:2])+'_xl.txt'
    target_dir =  '../../04_topo/'+symm_combo_flip+'/'
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for i,pdb in enumerate(pdb_list):
        print ('='*50)
        os.chdir(pdb_dir)
        if i == 0:
           write_xl(n,pdb_dir+pdb,'_'.join(symm_combo_flip.split('_')[:2]),xl_dir)
        Jwalk_run(jwalk_path,pdb,xl_dir,xl_txt,target_dir)
        shutil.move('./Jwalk_results/'+pdb[:-4]+'_crosslinks.pdb',target_dir+pdb[:-4]+xl_txt[9:-12]+'_crosslinks.pdb')
        shutil.move('./Jwalk_results/'+pdb[:-4]+'_crosslink_list.txt',target_dir+pdb[:-4]+xl_txt[9:-12]+'_crosslink_list.txt')
        
        xl_result_txt = target_dir+pdb[:-4]+xl_txt[9:-12]+'_crosslink_list.txt'
        xl_result_pdb = target_dir+pdb[:-4]+xl_txt[9:-12]+'_crosslinks.pdb'
        xl_temp_pdb = target_dir+'xl_tmp.pdb'
        altered_pdb = pdb[:-4]+xl_txt[9:-12]+'_crosslink_combined_semi.pdb'
        altered_pdb_f =  target_dir+altered_pdb
        os.chdir(target_dir)
        parse_crosslink(n,pdb[:-4]+xl_txt[9:-12],xl_result_txt,xl_result_pdb,xl_temp_pdb,pdb_dir+pdb,altered_pdb_f,xl_dir,target_dir)
        shutil.move('./Jwalk_results/'+altered_pdb[:-4]+'_crosslinks.pdb',target_dir+altered_pdb[:-4]+xl_txt[9:-12]+'_crosslinks.pdb')
        shutil.move('./Jwalk_results/'+altered_pdb[:-4]+'_crosslink_list.txt',target_dir+altered_pdb[:-4]+xl_txt[9:-12]+'_crosslink_list.txt')

        os.chdir(target_dir)
        xl_second_result_txt = target_dir+altered_pdb[:-4]+xl_txt[9:-12]+'_crosslink_list.txt'
        xl_second_result_pdb = target_dir+altered_pdb[:-4]+xl_txt[9:-12]+'_crosslinks.pdb'
        xl_second_temp_pdb = target_dir+'xl_second_tmp.pdb'
        combo_semi_pdb = altered_pdb_f
        altered_semi_pdb = target_dir+pdb[:-4]+xl_txt[9:-12]+'_crosslink_combined_whole.pdb'
        parse_second_crosslink(n,xl_second_result_txt,xl_second_result_pdb,xl_second_temp_pdb,combo_semi_pdb,altered_semi_pdb)