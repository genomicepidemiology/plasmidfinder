PlasmidFinder
===================

This project documents the PlasmidFinder service


Documentation
=============

The PlasmidFinder service contains one python script *plasmidfinder.py* which is the script of the latest
version of the PlasmidFinder service. The service identifies plasmids in total or partial sequenced
isolates of bacteria.

## Content of the repository
1. plasmidfinder.py      - the program
2. test     	- test folder
3. README.md
4. Dockerfile   - dockerfile for building the plasmidfinder docker container


## Installation

Setting up PlasmidFinder program
```bash
# Go to wanted location for plasmidfinder
cd /path/to/some/dir
# Clone and enter the plasmidfinder directory
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
cd plasmidfinder
```

Build Docker container
```bash
# Build container
docker build -t plasmidfinder .
# Run test
docker run --rm -it \
       --entrypoint=/test/test.sh plasmidfinder
```

#Download and install PlasmidFinder database

```bash
# Go to the directory where you want to store the plasmidfinder database
cd /path/to/some/dir
# Clone database from git repository (develop branch)
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
cd plasmidfinder_db
PLASMID_DB=$(pwd)
# Install PlasmidFinder database with executable kma_index program
python3 INSTALL.py kma_index
```

If kma_index has no bin install please install kma_index from the kma repository:
https://bitbucket.org/genomicepidemiology/kma

## Dependencies
In order to run the program without using docker, Python 3.5 (or newer) should be installed along with the following versions of the modules (or newer).

#### Modules
- cgecore 1.5.5
- tabulate 0.7.7

Modules can be installed using the following command. Here, the installation of the module cgecore is used as an example:
```bash
pip3 install cgecore
```
#### KMA and BLAST
Additionally KMA and BLAST version 2.8.1 or newer should be installed.
The newest version of KMA and BLAST can be installed from here:
```url
https://bitbucket.org/genomicepidemiology/kma
```

```url
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

## Usage

The program can be invoked with the -h option to get help and more information of the service.
Run Docker container:


```bash
# Run plasmidfinder container
docker run --rm -it \
       -v $PLASMID_DB:/database \
       -v $(pwd):/workdir \
       plasmidfinder -i [INPUTFILE] -o [OUTDIR] [-d] [-p] [-mp] [-l] [-t]
```

When running the docker file you have to mount 2 directories: 
 1. plasmidfinder_db (PlasmidFinder database) downloaded from bitbucket
 2. An output/input folder from where the input file can be reached and an output files can be saved. 
Here we mount the current working directory (using $pwd) and use this as the output directory, 
the input file should be reachable from this directory as well. The path to the infile and outfile
directories should be relative to the monuted current working directory.


`-i INPUTFILE	input file (fasta or fastq) relative to pwd, up to 2 files`

`-o OUTDIR	outpur directory relative to pwd`

`-d DATABASE    set a specific database`

`-p DATABASE_PATH    set path to database, default is /database`

`-mp METHOD_PATH    set path to method (blast or kma)`

`-l MIN_COV    set threshold for minimum coverage`

`-t THRESHOLD set threshold for mininum blast identity`


## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/PlasmidFinder/

Citation
=======

When using the method please cite:

PlasmidFinder and pMLST: in silico detection and typing of plasmids.
Carattoli A, Zankari E, Garcia-Fernandez A, Volby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.
Antimicrob. Agents Chemother. 2014. April 28th.
[Epub ahead of print]

References
=======

1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics 2009; 10:421. 
2. Clausen PTLC, Aarestrup FM, Lund O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 2018; 19:307. 

License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.