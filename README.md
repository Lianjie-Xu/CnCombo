# CnCombo
This is the code for protein topology design from the combination of Cn entangling motifs as described in the paper "*Computational Design and Cellular Synthesis of Two Protein Topological Isoforms: Solomon Link vs Three-Twist Knot*".

# Installation
Make a conda environment to run *CnCombo* (PyRosetta license is required to fill in USERNAME:PASSWORD in environment.yml):
```shell
git clone https://github.com/Lianjie-Xu/CnCombo.git     
cd CnCombo    
conda env create --name=CnCombo -f environment.yml    
conda activate CnCombo
```
-----------------------------------------------------------------------------------------------------
Code organization:
* `CnCombo_run.py` - the main script to run simulation.
* `combo_multimer.py` - script to prepare Cn entangling motifs as inputs.
* `combo_symm.py` - script to generate conformations of a given combination of Cn motifs.
* `combo_trans_rot_Cn.py` - script to sample conformations.
* `combo_linkages_Cn.py` - script to calculate minimum SASDs for the conformational ensemble.
* `combo_topo_Cn.py` - script to determine topology.
* `helper_scripts/` - helper softwares to handle symmetry and to calculate solvent-accessible surface distances.
* `motifs/` - collection of Cn entangling motifs.
* `outputs/` - outputs from simulation.
-----------------------------------------------------------------------------------------------------
Input flags for `CnCombo_run.py`:
```
    argparser.add_argument("--suppress_print", type=int, default=0, help="0 for False, 1 for True")
    parser.add_argument('-i1', '--input_pdb1', type=str, help='input pdb file of motif1')
    parser.add_argument('-i2', '--input_pdb2', type=str, help='input pdb file of motif2')
    parser.add_argument('-t', '--target_axis', type=str, default='0,0,1', help='axis to be aligned to')
    parser.add_argument('-symm1', '--symm1', type=str, help='symmetry symbol of input motif1')
    parser.add_argument('-symm2', '--symm2', type=str, help='symmetry symbol of input motif2')
    parser.add_argument('-s', '--target_symm', type=str, help='target assembly symmetry')
    parser.add_argument('-c', '--cutoff', type=float, default=5.0, help='distance cutoff for motif translation')
    parser.add_argument('-g', '--interval', type=float, default=5.0, help='interval for motif rotation')
```
-----------------------------------------------------------------------------------------------------
# Usage (taking X+K as an example)   
```shell
python CnCombo_run.py --input_pdb1 ./motifs/1sak.pdb --symm1 C2 --input_pdb2 ./motifs/2bo3.pdb --symm2 C2 --target_symm C2
```

-----------------------------------------------------------------------------------------------------
```
@article{
  title={Computational Design and Cellular Synthesis of Two Protein Topological Isoforms: Solomon Link vs Three-Twist Knot},
  author={xxxx},
  journal={xxxx},
  volume={xxxx},
  number={xxxx},  
  pages={xxx--xx},
  year={xxxx},
  publisher={xxxxxx}
}
```
-----------------------------------------------------------------------------------------------------
