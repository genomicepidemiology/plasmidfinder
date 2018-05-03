#!/usr/bin/env python2.7
from __future__ import division
import sys, os, time, random, re, subprocess
from argparse import ArgumentParser
from tabulate import tabulate
import collections
from blaster import *
from distutils.spawn import find_executable

##########################################################################
# FUNCTIONS
##########################################################################

def text_table(headers, rows, empty_replace='-'):
   ''' Create text table
   
   USAGE:
      >>> from tabulate import tabulate
      >>> headers = ['A','B']
      >>> rows = [[1,2],[3,4]]
      >>> print(text_table(headers, rows))
      **********
        A    B
      **********
        1    2
        3    4
      ==========
   '''
   # Replace empty cells with placeholder
   rows = map(lambda row: map(lambda x: x if x else empty_replace, row), rows)
   # Create table
   table = tabulate(rows, headers, tablefmt='simple').split('\n')
   # Prepare title injection
   width = len(table[0])
   # Switch horisontal line
   table[1] = '*'*(width+2)
   # Update table with title
   table = ("%s\n"*3)%('*'*(width+2), '\n'.join(table), '='*(width+2))
   return table

##########################################################################
#	DEFINE GLOBAL VARIABLES
##########################################################################

global database_path, databases, min_length, threshold

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", dest="inputfile",help="Input file", default='')
parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output", default='')
parser.add_argument("-b", "--blastPath", dest="blast_path",help="Path to blast", default='blastn')
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='')
parser.add_argument("-d", "--databases", dest="databases",help="Databases chosen to search in - if non is specified all is used", default=None)
parser.add_argument("-l", "--mincov", dest="min_cov",help="Minimum coverage", default=0.60)
parser.add_argument("-t", "--threshold", dest="threshold",help="Blast threshold for identity", default=0.90)
args = parser.parse_args()


##########################################################################
# MAIN
##########################################################################

# Defining varibales

min_cov = args.min_cov
threshold = args.threshold

# Check if valid database is provided
if args.db_path is None:
      sys.exit("Input Error: No database directory was provided!\n")
elif not os.path.exists(args.db_path):
   sys.exit("Input Error: The specified database directory does not"
                       " exist!\n")
else:
   # Check existence of config file
   db_config_file = '%s/config'%(args.db_path)
   if not os.path.exists(db_config_file):
      sys.exit("Input Error: The database config file could not be "
                          "found!")
   # Save path
   db_path = args.db_path

# Check if valid input file is provided
if args.inputfile is None:
   sys.exit("Input Error: No Input were provided!\n")
elif not os.path.exists(args.inputfile):
   sys.exit("Input Error: Input file does not exist!\n")
else:
    inputfile = args.inputfile

# Check if valid output directory is provided
if not os.path.exists(args.out_path):
   # sys.exit("Input Error: Output dirctory does not exists!\n")
   out_path = '.'
else:
   out_path = args.out_path

# Check if valid path to BLAST is provided
if not os.path.exists(args.blast_path):
   blast = 'blastn'
   if find_executable(blast) is None:
      sys.exit("Input Error: blastn could not be found, please provide a path to blastn!\n")
else:
   blast = args.blast_path

# Check if databases and config file are correct/correponds
if args.databases is '':
      sys.exit("Input Error: No database was specified!\n")
else:
   dbs = dict()
   extensions = []
   with open(db_config_file) as f:
      for l in f:
         l = l.strip()
         if l == '': continue
         if l[0] == '#':
            if 'extensions:' in l:
               extensions = [s.strip() for s in l.split('extensions:')[-1].split(',')]
            continue
         tmp = l.split('\t')
         if len(tmp) != 3:
            sys.exit(("Input Error: Invalid line in the database"
                                 " config file!\nA proper entry requires 3 tab "
                                 "separated columns!\n%s")%(l))
         db_prefix = tmp[0].strip()
         name = tmp[1].split('#')[0].strip()
         description = tmp[2]
         # Check if all db files are present
         for ext in extensions:
            db = "%s/%s.%s"%(db_path, db_prefix, ext)
            if not os.path.exists(db):
               sys.exit(("Input Error: The database file (%s) "
                                    "could not be found!")%(db_path))
         if db_prefix not in dbs: dbs[db_prefix] = []
         dbs[db_prefix].append(name)
   if len(dbs) == 0:
      sys.exit("Input Error: No databases were found in the "
                          "database config file!")

   if args.databases is None:
      # Choose all available databases from the config file
      databases = dbs.keys()
   else:
      # Handle multiple databases
      args.databases = args.databases.split(',')
      # Check that the ResFinder DBs are valid
      databases = []
      for db_prefix in args.databases:
         if db_prefix in dbs:
            databases.append(db_prefix)
         else:
            sys.exit("Input Error: Provided database was not "
                     "recognised! (%s)\n"%db_prefix)

# Calling blast and parsing output
results, query_align, homo_align, sbjct_align = Blaster(
   inputfile, databases, db_path, out_path, min_cov, threshold, blast)

# Getting and writing out the results
rows = list()
dbs_with_results = list()
txt_file_seq_text = dict()
split_list = collections.defaultdict(list)
join_as_str = lambda x, y='\t': y.join(map(str, x))

headers = ["Database", "Plasmid", "Identity", "Alignment Length", "Template Length", "Position in reference", "Contig", "Position in contig", "Note", "Accession no."]
with open(out_path+"/results_tab.txt", 'w') as tab_file, \
     open(out_path+"/results_table.txt", 'w') as table_file, \
     open(out_path+"/Plasmid_seq.fsa", 'w') as ref_file, \
     open(out_path+"/Hit_in_genome_seq.fsa", 'w') as hit_file:
   # Write the header for the tab file
   tab_file.write("%s\n"%(join_as_str(headers)))
   table_file.write("%s\n"%(join_as_str(headers)))
   
   for db in results:
      db_name = str(dbs[db][0])
      
      if results[db] == "No hit found":
         table_file.write("%s\nNo hits found\n\n"%db_name)
      else:
         dbs_with_results.append(db_name)
         txt_file_seq_text[db_name] = list()
         for hit in results[db]:
            header = results[db][hit]["sbjct_header"]
            tmp = header.split("_")
            gene = tmp[0]
            note = tmp[2]
            acc = tmp[3]
            ID = results[db][hit]["perc_ident"]
            sbjt_length = results[db][hit]["sbjct_length"]
            HSP = results[db][hit]["HSP_length"]
            positions_contig = "%s..%s"%(results[db][hit]["query_start"], results[db][hit]["query_end"])
            positions_ref = "%s..%s"%(results[db][hit]["sbjct_start"], results[db][hit]["sbjct_end"])
            contig_name = results[db][hit]["contig_name"]
            
            tmp_arr = [db_name, gene, ID, HSP, sbjt_length, positions_ref, contig_name,
                     positions_contig, note, acc]
            if "split_length" in results[db][hit]:
               tab_file.write("%s\n"%(join_as_str(tmp_arr)))
               
               # Switch hsp_length with total HSP and push to split list
               tmp_arr[3] = results[db][hit]["split_length"]
               split_list[res].append(tmp_arr)
            else:
               # Saving the output to write the txt result table
               rows.append(tmp_arr[:])
               
               # Write tabels
               tmp_arr[2] = "%.2f"%tmp_arr[2]
               tmp = "%s\n"%(join_as_str(tmp_arr))
               table_file.write(tmp)
               tab_file.write(tmp)
            
            # Writing subjet/ref sequence
            ref_seq = sbjct_align[db][hit]
            ref_file.write(">%s_%s\n"%(gene, acc))
            for i in range(0, len(ref_seq), 60):
               ref_file.write("%s\n"%(ref_seq[i:i + 60]))
            
            # Getting the header and text for the txt file output
            #>aac(2')-Ic: PERFECT MATCH, ID: 100.00%, HSP/Length: 546/546, Positions in reference: 1..546, Contig name: gi|375294201|ref|NC_016768.1|, Position: 314249..314794
            sbjct_start = results[db][hit]["sbjct_start"]
            sbjct_end = results[db][hit]["sbjct_end"]
            text = "%s, ID: %.2f %%, Alignment Length: %s Template Length: %s, Positions in reference: %s..%s, Contig name: %s, Position: %s"%(gene, ID, HSP, sbjt_length, sbjct_start, sbjct_end, contig_name, positions_contig)
            hit_file.write(">%s\n"%text)
            
            # Writing query/hit sequence
            hit_seq = query_align[db][hit]
            for i in range(0, len(hit_seq), 60):
               hit_file.write("%s\n"%(hit_seq[i:i + 60]))
            
            # Saving the output to print the txt result file allignemts
            txt_file_seq_text[db_name].append((text, ref_seq, homo_align[db][hit], hit_seq))
         
         for res in split_list:
            gene = split_list[res][0][0]
            ID = split_list[res][0][1]
            HSP = split_list[res][0][2]
            sbjt_length = split_list[res][0][3]
            positions_ref = split_list[res][0][4]
            contig_name = split_list[res][0][5]
            positions_contig = split_list[res][0][6]
            note = split_list[res][0][7]
            acc = split_list[res][0][8]
            
            for i in range(1,len(split_list[res])):
               ID = "%s, %.2f"%(ID, split_list[res][i][1])
               positions_ref = positions_ref + ", " + split_list[res][i][4]
               contig_name = contig_name + ", " + split_list[res][i][5]
               positions_contig = positions_contig + ", " + split_list[res][i][6]
            
            tmp_arr = [db_name, gene, ID, HSP, sbjt_length,
                       positions_ref, contig_name, positions_contig, note, acc]
            table_file.write("%s\n"%(join_as_str(tmp_arr)))
            print(tmp_arr)
            rows.append(tmp_arr)
         
         table_file.write("\n")

# Writing the txt file
with open(out_path+"/results.txt", 'w') as f:
   # Writing table
   f.write(text_table(headers, rows))
   
   # Writing extended output
   if dbs_with_results:
      f.write("\n\nExtended Output:\n")
      f.write("%s\n"%('-'*80))
      for db_name in dbs_with_results:
         # Txt file alignments
         # test = lambda x: (x, (80 - x)//2, x%2, ((80 - x)//2 - 2) * 2  + (1 if x%2 else 0) + x +4)
         tlen = len(db_name)
         spacer_size = (80 - tlen)//2 - 2
         spacer_right = '#'*spacer_size
         if tlen%2: spacer_size += 1
         
         spacer_left = '#'*spacer_size
         f.write("#%s %s %s#\n"%(spacer_left, db_name, spacer_right))
         for text in txt_file_seq_text[db_name]:
            f.write("%s\n\n"%(text[0]))
            for i in range(0, len(text[1]), 60):
               f.write("%-20s%s\n"%('Template:', text[1][i:i + 60]))
               f.write("%-20s%s\n"%('', text[2][i:i + 60]))
               f.write("%-20s%s\n\n"%('Query:', text[3][i:i + 60]))
            
            f.write("\n%s\n"%('-'*80))
