# Endothermy project

This repo consists of the scripts used for our project studying the evolution of endothermy. We outline the selection analyses, for other analyses see dryad. 

This project uses various methods from the [HyPhy package](https://github.com/veg/hyphy). For specific batch files (.bf), see HyPhy analyses repo. 

## The pipeline
Our pipeline was developed to be used in an HPC environment. 
There is an assumption that the freely availible [Anaconda](https://anaconda.org/) software is installed on your machine. 

### To install:

1. `git clone https://github.com/agselberg/Endothermy.git`
2. `cd Endothermy`
3. `conda env create -f environment.yml`.  This will create a conda environment with the necessary dependencies.
4. At this point, run `conda activate ENDOTHERMY` and your environment will be ready to go.
