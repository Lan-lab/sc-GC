# Single cell analysis coupled with TCR sequencing of gastric tumor reveals complex intercellular interaction and an alternative T cell exhaustion trajectory
Code for the central analyses in the study
## Environment 
Red Hat 4.8.5-16, Linux version 3.10.0
R version 3.5.1	
Python version 3.6.7
seaborn-0.11.1
## Install software
This should take less than an hour, with fast internet speed (e.g. 100 Mbps).

### Install basic python packages 
```
pip install numpy==1.19.5 pandas==1.1.5 matplotlib==3.3.4 seaborn==0.11.1 scipy==1.5.4
```
### Install python package scanpy 
```
pip install scanpy==1.7.2
```
### Install python package scvelo 
```
pip install scvelo=0.2.3
```
### Install python package pySCENIC
    pip install pyscenic==0.11.0
### Install R package InferCNV 
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("infercnv")
### Insatll software CellPhoneDB v2.06
    pip install cellphonedb==2.06
	

## Possible issue

### If there is any error with the umap module inside scanpy, reinstall umap-learn with conda
```
pip uninstall umap-learn
conda install umap-learn
```


## Runing the code
	The scripts with .py and .r should be tested in an interactive console.
	The run time for demo in 'Data' should be several minutes for each analysis.