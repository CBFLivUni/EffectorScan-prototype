# EffectorScan-prototype
 
Prototype code for identifying novel effectors from fusarium oxysporum unique genes. It starts from the raw fasta files 'core_Fo_genes_with_SP.fasta' and 'lycopersici_unique_with_SP.fasta' contained in a directory named 'raw_fasta'. The code is a mixture of bash scripts and Quarto notebooks (with R script alternatives available for the Quarto notebooks). These should be run in order and contain documentation within. 

It is reccomended this workflow is ran on a cluster with a GPU to speed up the computational and creation of structure models using colabfold. 

The analysis requires two key pieces of software to be installed (in addition to a standard Rstudio set-up): 
- [ColabFold](https://github.com/sokrypton/ColabFold) installation instructions can be found here
- [FoldSeek](https://github.com/steineggerlab/foldseek) installation instructions and tutorials can be found here