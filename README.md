# SolCuration

This is the repo for Tencent SolCuration Project. We provide an Aqueous Solubility Dataset Curation tool and 7 Solubility Dataset Before and After Data Curation.


## HOWTO

### How to download SolCuration tools and datasets  

Simply run `git clone https://github.com/Mengjintao/SolCuration`.

The above command will download the data curation code and all the original and curated datasets automatically.  You can also download your favorite dataset directly from our repo manually.  

### How to evaluate the dataset with Chemprop and AttentiveFP

#### Setup the Conda virtual Enviroment

The easiest way to install the dependencies for Chemprop and AttentiveFP is via conda. Here are the steps:

1. Install conda from [Anaconda Distribution](https://docs.anaconda.com/anaconda/install/)
2. `cd /path/to/SolCuration`
3. `conda env create -f environment.yml`
4. `source activate mpnn` (or `conda activate mpnn` for newer versions of conda)
5. `pip install -r requirements.txt`

After the installation, you can check your enviroment with the following commands:

1. `conda activate mpnn`
2. `python --version`  should returns `Python 3.7.6`
3. `python` 
4. `import torch`
5. `print (torch.cuda.device_count())`  should returns the number of GPU cards in your host. 

#### Evaluate with Chemprop

1. `cd /path/to/SolCuration`
2. `git clone https://github.com/Mengjintao/Chemprop`
3. `cd Chemprop`
4. `bash test.sh` can lanch a test run for Chemprop.

#### Evaluate with AttentiveFP

1. `cd /path/to/SolCuration`
2. `git clone https://github.com/Mengjintao/AttentiveFP`
3. `cd Chemprop`
4. `python runScaffold.py ../../org/esol/esol_org.csv 124 0 2 5 100  5 2` can lanch a test run for AttentiveFP.

#### Running Hyper-Parameter Using Grid Search With Slurm On A GPU Cluster
Currently, we use szsc partition (2000 computing nodes with 6,400 CPU cores and 8000 GPUs) in National Supercomputer Center in Shenzhen as our GPU cluster. 

1. `cd /path/to/Chemprop`
2. `bash Chemprop_GridSearch.sh` 
3. `cd /path/to/AttentiveFP`
4. `bash AttentiveFP_GridSearch.sh` 

##### Collect RMSE and confidence intervals for 7 datasets After Grid Search
Chemprop
1. `cd /path/to/Chemprop`
2. `bash Collect_RMSE.sh > RMSE.txt`
3.  We can load data from RMSE.txt using MS Excel for further analysis.  

AttentiveFP
1. `cd /path/to/AttentiveFP`
2. `bash Collect_RMSE.sh > RMSE.txt`
3.  We can load data from RMSE.txt using MS Excel for further analysis. 

##### Collect Pearson Correlation Coefficient R^2 and Spearman's Rank-Order Correlation Coefficient R_s on BPU, BDZ, PCA, CDK, BPZ&BDZ series evaluation dataset.  
Chemprop
1. `cd /path/to/Chemprop`
2. `bash Collect_Coefficient.sh > R2out.txt`
3.  We can load data from R2out.txt using MS Excel and plot the figures for further analysis. 

AttentiveFP
1. `cd /path/to/AttentiveFP`
2. `bash Collect_Coefficient.sh > R2out.txt`
3.  We can load data from R2out.txt using MS Excel and plot the figures for further analysis. 

#### Dependencies  

+ Python 3.6+
+ PyTorch 1.0+
+ RDKit
+ torchvision
+ pandas
+ tqdm
+ openbabel

### HOWTO for Advanced Users

#### Adjust the data curation parameters and rerun the data curation workflow to generate customed curated datasets.

1. `cd /path/to/SolCuration/src`
2. Modify or add the weight at line 477~484 in cure.cpp.
3. recompile cure.cpp with `g++ cure.cpp -o cure`
4. `cd ..`
5. Generated customed datasets with `bash dataGeneration.sh`
6. Generated additional features with rdkit `bash featureGeneration.sh`

#### Reference  

- Sorkun, Murat Cihan and Khetan, Abhishek and Er, Suleyman, AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds, Nature Scientific data, 6, pp.1-8 (2019) [link](https://www.nature.com/articles/s41597-019-0151-1)
- Yang, Kevin and Swanson, Kyle and Jin, Wengong and Coley, Connor and Eiden, Philipp and Gao, Hua and Guzman-Perez, Angel and Hopper, Timothy and Kelley, Brian and Mathea, Miriam and others.
Analyzing learned molecular representations for property prediction. Journal of chemical information and modeling, 59, pp. 3370-3388 (2019) [link](https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.9b00237)  
- Xiong, Zhaoping and Wang, Dingyan and Liu, Xiaohong and Zhong, Feisheng and Wan, Xiaozhe and Li, Xutong and Li, Zhaojun and Luo, Xiaomin and Chen, Kaixian and Jiang, Hualiang and others, Pushing the boundaries of molecular representation for drug discovery with the graph attention mechanism. Journal of medicinal chemistry (2019) [link](https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b00959) 
 

## Citation

Please cite as:
```
```
