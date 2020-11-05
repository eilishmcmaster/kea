# Kea: MAG scaffolding and gap closing using draft genomes

The Kea scaffolding tool uses multiple fragmented MAGs of the same OTU to scaffold contigs together, closing gaps in a representative genome. It is made to be used in a Linux command line.


## Installation 
The Kea program is still in beta development and hasn't been wrapped yet. You will need to git clone and create a conda environment before use. 

### Install Kea scripts

1. clone this git repository. This will create the directory `kea` containing all scripts.

```
git clone https://github.com/eilishmcmaster/kea.git
```

2. Activate the python scripts
```
chmod +x /kea/kea_modular/*py
```

### Create conda environment 

This environment will contain all of the programs needed for Kea to run 

1. Make new environment 

```
conda create -y -p kea-env
```

2. Activate the environment (command should be dispayed at the end of the `conda create` output)

```
conda activate kea-env
```

3. Install the prerequisite packages into the environment

```
conda install --name kea-env --file /kea/conda_install.txt
```
All the prerequisite packages should now be installed. 

## Using Kea 

To use Kea, you must have multiple MAGs of the same OTU (>99% ANI). If you have multiple different OTUs you would like to run Kea on, Kea can sort them automatically. The input is the directory containing the MAGs. 

Input format: 

```
python kea/kea_modular/main.py -i input_directory/ -o output_directory -x fa -r reference.fa -t 50
```

The `x` flag is used to specify the input file extension (fa or fasta), `t` is to specify the number of threads to use, and `r` is for specifying the reference MAG (optional). The output file will be created if it does not exist. 
