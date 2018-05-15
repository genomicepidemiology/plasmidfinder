===================
PlasmidFinder
===================

This project documents PlasmidFinder service


Documentation
=============

## What is it?

The PlasmidFinder service contains one perl script *plasmidfinder.pl* which is the script of the lates
version of the PlasmidFinder service. The service identifies plasmids in total or partial sequenced
isolates of bacteria.

## Content of the repository
1. plasmidfinder.pl  - the program
2. INSTALL_DB        - shell script for downloading the MLST database
3. UPDATE_DB         - shell script for updating the database to the newest version
4. VALIDATE_DB       - python script for verifying the database contains all
                       required files 
5. brew.sh           - shell script for installing dependencies
6. makefile          - make script for installing dependencies
7. test.fsa          - test fasta file

## Installation

### Setting up PlasmidFinder
```bash
# Go to wanted location for resfinder
cd /path/to/some/dir
# Clone and enter the mlst directory
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
cd plasmidfinder
```

### Installing up the PlasmidFinder database
Check out the database repository.
[https://bitbucket.org/genomicepidemiology/plasmidfinder_db]

### Installing dependencies:

Perlbrew is used to manage isolated perl environments. To install it run:
```bash
bash brew.sh
```

This will installed Perl 5.23 in the Home folder, along with CPAN minus as package manager.
Blast will also be installed when running brew.sh if BlastAll and FormatDB are not already installed and place in the user's path.
After running brew.sh and installing Blast add this command to the end of your ~/bash_profile to add BlastAll and FormatDB to the user's path

```bash
export PATH=$PATH:blast-2.2.26/bin
```

If you want to download the two external tools from the Blast package, BlastAll and FormatDB, yourself go to
```url
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST
```

and download the version for your OS with the format:
```url
blast-version-architecture-OS.tar.gz
```

after unzipping the file, add this command to the end of your ~/bash_profile.
```bash
export PATH=$PATH:/path/to/blast-folder/bin
```

where path/to/blast-folder is the folder you unzipped.

At last PlasmidFinder has several Perl dependencies. To install them (this requires CPAN minus as package manager):
```bash
make install
```

The scripts are self contained. You just have to copy them to where they should
be used. Only the *database* folder needs to be updated mannually.

Remember to add the program to your system path if you want to be able to invoke the program without calling the full path.
If you don't do that you have to write the full path to the program when using it.

## Usage

The program can be invoked with the -h option to get help and more information of the service.

```bash
Usage: perl plasmidfinder.pl [options]

Options:

    -h HELP
                    Prints a message with options and information to the screen
    -d DATABASE
                    The path to where you have located the database folder
    -b BLAST
                    The path to the location of blast-2.2.26 if it is not added
                    to the users path (see the install guide in 'README.md')
    -i INFILE
                    Your input file which needs to be preassembled partial
                    or complete genomes in fasta format
    -o OUTFOLDER
                    The folder you want to have your output files places.
                    If not specified the program will create a folder named
                    'Output' in which the result files will be stored.
    -p PLASMID_DATBASE
                    The database you want to search for the plasmids in.
                    Choose ehiter plasmid_database, plasmid_positiv
                    or both (plasmid_database,plasmid_positiv)
    -k  THRESHOLD
                    The threshold for % identity for example '95.00' for 95 %
```

#### Example of use with the *database* folder located in the current directory and Blast added to the user's path
```perl
    perl plasmidfinder.pl -i test.fsa -o OUTFOLDER -p plasmid_database -k 95.00
```
#### Example of use with the *database* and *blast-2.2.26* folders loacted in other directories
```perl
    perl plasmidfinder.pl -d path/to/database -b path/to/blast-2.2.26 -i \
    test.fsa -o OUTFOLDER -p plasmid_database -k 95.00
```

## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/PlasmidFinder/


## The Latest Version


The latest version can be found at
https://bitbucket.org/genomicepidemiology/plasmidfinder/overview

## Documentation


The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/plasmidfinder/overview.

Citation
=======

When using the method please cite:

PlasmidFinder and pMLST: in silico detection and typing of plasmids.
Carattoli A, Zankari E, Garcia-Fernandez A, Volby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.
Antimicrob. Agents Chemother. 2014. April 28th.
[Epub ahead of print]


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
