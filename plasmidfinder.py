#!/usr/bin/env python3
from __future__ import division
import sys, os, time, random, re, subprocess
from argparse import ArgumentParser
from tabulate import tabulate
import collections
from cgecore.blaster import Blaster
from cgecore.cgefinder import CGEFinder
from distutils.spawn import find_executable
import json, gzip, pprint

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

def get_read_filename(infiles):
   ''' Infiles must be a list with 1 or 2 input files.
       Removes path from given string and removes extensions:
       .fq .fastq .gz and .trim
       extract the common sample name i 2 files are given.
   '''
   # Remove common fastq extensions
   seq_path = infiles[0]
   seq_file = os.path.basename(seq_path)
   seq_file = seq_file.replace(".fq", "")
   seq_file = seq_file.replace(".fastq", "")
   seq_file = seq_file.replace(".gz", "")
   seq_file = seq_file.replace(".trim", "")
   if len(infiles) == 1:
      return seq_file.rstrip()

   # If two files are given get the common sample name
   sample_name = ""
   seq_file_2 = os.path.basename(infiles[1])
   for i in range(len(seq_file)):
      if seq_file_2[i] == seq_file[i]:
         sample_name += seq_file[i]
      else: 
         break
   if sample_name == "":
      sys.error("Input error: sample names of input files, {} and {}, \
                 does not share a common sample name. If these files \
                 are paired end reads from the same sample, please rename \
                 them with a common sample name (e.g. 's22_R1.fq', 's22_R2.fq') \
                 or input them seperately.".format(infiles[0], infiles[1]))

   return sample_name.rstrip("-").rstrip("_")

def is_gzipped(file_path):
   ''' Returns True if file is gzipped and False otherwise.
       The result is inferred from the first two bits in the file read
       from the input path.
       On unix systems this should be: 1f 8b
       Theoretically there could be exceptions to this test but it is
       unlikely and impossible if the input files are otherwise expected
       to be encoded in utf-8.
   '''
   with open(file_path, mode='rb') as fh:
      bit_start = fh.read(2)
   if(bit_start == b'\x1f\x8b'):
      return True
   else:
      return False

def get_file_format(input_files):
   """
   Takes all input files and checks their first character to assess
   the file format. Returns one of the following strings; fasta, fastq, 
   other or mixed. fasta and fastq indicates that all input files are 
   of the same format, either fasta or fastq. other indiates that all
   files are not fasta nor fastq files. mixed indicates that the inputfiles
   are a mix of different file formats.
   """

   # Open all input files and get the first character
   file_format = []
   invalid_files = []
   for infile in input_files:
      if is_gzipped(infile):#[-3:] == ".gz":
         f = gzip.open(infile, "rb")
         fst_char = f.read(1);
      else:
         f = open(infile, "rb")
         fst_char = f.read(1);
      f.close()
      # Assess the first character
      if fst_char == b"@":
         file_format.append("fastq")
      elif fst_char == b">":
         file_format.append("fasta")
      else:
         invalid_files.append("other")
   if len(set(file_format)) != 1:
      return "mixed"
   return ",".join(set(file_format))

##########################################################################
#	DEFINE GLOBAL VARIABLES
##########################################################################

global database_path, databases, min_length, threshold

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", dest="inputfile",help="Input file", default='', nargs="+")
parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output", default='')
parser.add_argument("-tmp", "--tmp_dir")
parser.add_argument("-mp", "--methodPath", dest="method_path",help="Path to method", default='')
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='/database')
parser.add_argument("-d", "--databases", dest="databases",help="Databases chosen to search in - if non is specified all is used", default=None)
parser.add_argument("-l", "--mincov", dest="min_cov",help="Minimum coverage", default=0.60)
parser.add_argument("-t", "--threshold", dest="threshold",help="Blast threshold for identity", default=0.90)
parser.add_argument("-x", "--extented_output",
                    help="Give extented output with allignment files, template and query hits in fasta and\
                          a tab seperated file with allele profile results", action="store_true")

args = parser.parse_args()


##########################################################################
# MAIN
##########################################################################

# Defining varibales

min_cov = float(args.min_cov)
threshold = float(args.threshold)
method_path = args.method_path

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

# Check if valid input files are provided
if args.inputfile is None:
   sys.exit("Input Error: No Input were provided!\n")

elif not os.path.exists(args.inputfile[0]):
   sys.exit("Input Error: Input file does not exist!\n")
elif len(args.inputfile) > 1 and not os.path.exists(args.inputfile[1]):
   sys.exit("Input Error: Input file does not exist!\n")

else:   
   inputfile = args.inputfile

# Check if valid output directory is provided
if not os.path.exists(args.out_path):
   # sys.exit("Input Error: Output dirctory does not exists!\n")
   out_path = '.'
else:
   out_path = os.path.abspath(args.out_path)

# Check if valid tmp directory is provided
if args.tmp_dir:
   if not os.path.exists(args.tmp_dir):
      sys.exit("Input Error: Tmp dirctory, {}, does not exists!\n".format(args.tmp_dir))
   else:
      tmp_dir = os.path.abspath(args.tmp_dir)
else:
   tmp_dir = out_path

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

# Check file format (fasta, fastq or other format)
file_format = get_file_format(inputfile)

# Call appropriate method (kma or blastn) based on file format 
if file_format == "fastq":
   if not method_path:
      method_path = "kma"
      if find_executable(method_path) == None:
         sys.exit("No valid path to a kma program was provided. Use the -mp flag to provide the path.")
   # Check the number of files
   if len(inputfile) == 1:
      infile_1 = inputfile[0]
      infile_2 = None
   elif len(inputfile) == 2:
      infile_1 = inputfile[0]
      infile_2 = inputfile[1]
   else:
      sys.exit("Only 2 input file accepted for raw read data,\
                if data from more runs is avaliable for the same\
                sample, please concatinate the reads into two files")
    
   sample_name = get_read_filename(inputfile)
   method = "kma"

   # Call KMA
   method_obj = CGEFinder.kma(infile_1, tmp_dir, databases, db_path, min_cov=min_cov,
                              threshold=threshold, kma_path=method_path, sample_name=sample_name,
                              inputfile_2=infile_2, kma_mrs=0.75, kma_gapopen=-5,
                              kma_gapextend=-1, kma_penalty=-3, kma_reward=1)

elif file_format == "fasta":
   if not method_path:
      method_path = "blastn"
      if find_executable(method_path) == None:
         sys.exit("No valid path to a blastn program was provided. Use the -mp flag to provide the path.")
   # Assert that only one fasta file is inputted
   assert len(inputfile) == 1, "Only one input file accepted for assembled data"
   inputfile = inputfile[0]
   method = "blast"

   # Call BLASTn
   method_obj = Blaster(inputfile, databases, db_path, out_path, min_cov, threshold, method_path)
else:
   sys.exit("Input file must be fastq or fasta format, not "+ file_format)

results     = method_obj.results
query_align = method_obj.gene_align_query
homo_align  = method_obj.gene_align_homo
sbjct_align = method_obj.gene_align_sbjct

# Getting and writing out the results
rows = list()
dbs_with_results = list()
txt_file_seq_text = dict()
split_list = collections.defaultdict(list)
join_as_str = lambda x, y='\t': y.join(map(str, x))
json_results = dict()

tab_file_lst = []
table_file_lst = []
ref_file_lst = []
hit_file_lst = []

headers = ["Database", "Plasmid", "Identity", "Alignment Length", "Template Length", "Position in reference", "Contig", "Position in contig", "Note", "Accession no."]
header_str = "\t".join(headers) + "\n"

# Write the header for the tab files
tab_file_lst.append(header_str)
table_file_lst.append(header_str)
   
for db in results:
   if db == 'excluded':
      continue
   db_name = str(dbs[db][0])
   if db_name not in json_results:
      json_results[db_name] = {}
   if results[db] == "No hit found":
      table_file_lst.append("%s\nNo hits found\n\n"%db_name)
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
            tab_file_lst.append("%s\n"%(join_as_str(tmp_arr)))
            
            # Switch hsp_length with total HSP and push to split list
            tmp_arr[3] = results[db][hit]["split_length"]
            split_list[hit].append(tmp_arr)
         else:
            # Saving the output to write the txt result table
            rows.append(tmp_arr[:])
            
            # Write tabels
            tmp_arr[2] = "%.2f"%tmp_arr[2]
            tmp = "%s\n"%(join_as_str(tmp_arr))
            table_file_lst.append(tmp)
            tab_file_lst.append(tmp)
            
         # Writing subjet/ref sequence
         ref_seq = sbjct_align[db][hit]
         ref_file_lst.append(">%s_%s\n"%(gene, acc))
         for i in range(0, len(ref_seq), 60):
            ref_file_lst.append("%s\n"%(ref_seq[i:i + 60]))
           
         # Getting the header and text for the txt file output
         #>aac(2')-Ic: PERFECT MATCH, ID: 100.00%, HSP/Length: 546/546, Positions in reference: 1..546, Contig name: gi|375294201|ref|NC_016768.1|, Position: 314249..314794
         sbjct_start = results[db][hit]["sbjct_start"]
         sbjct_end = results[db][hit]["sbjct_end"]
         text = "%s, ID: %.2f %%, Alignment Length: %s, Template Length: %s, Positions in reference: %s..%s, Contig name: %s, Position: %s"%(gene, ID, HSP, sbjt_length, sbjct_start, sbjct_end, contig_name, positions_contig)

         hit_file_lst.append(">%s\n"%text)
         # Writing query/hit sequence
         hit_seq = query_align[db][hit]
         for i in range(0, len(hit_seq), 60):
            hit_file_lst.append("%s\n"%(hit_seq[i:i + 60]))
            
         # Saving the output to print the txt result file allignemts
         txt_file_seq_text[db_name].append((text, ref_seq, homo_align[db][hit], hit_seq))
      
         # Write JSON results dict
         json_results[db_name].update({header:{}})
         json_results[db_name][header] = {"plasmid":gene,"ID":round(ID,2),"align_length":HSP,
                                             "template_length":sbjt_length,"position_in_ref":positions_ref,
                                             "contig_name":contig_name,"positions_in_contig":positions_contig,
                                             "note":note,"accession":acc}
      # TODO What are the split genes? Are they saved in the json output?
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

         table_file_lst.append("%s\n"%(join_as_str(tmp_arr)))
         rows.append(tmp_arr)

      table_file_lst.append("\n")


# Get run info for JSON file
service = os.path.basename(__file__).replace(".py", "")
date = time.strftime("%d.%m.%Y")
time = time.strftime("%H:%M:%S")

# Make JSON output file
data = {service:{}}

userinput = {"filename(s)":args.inputfile, "method":method,"file_format":file_format}
run_info = {"date":date, "time":time,"dbs_with_results":dbs_with_results}

data[service]["user_input"] = userinput
data[service]["run_info"] = run_info
data[service]["results"] = json_results

pprint.pprint(data)

# Save json output
result_file = "{}/data.json".format(out_path) 
with open(result_file, "w") as outfile:
   json.dump(data, outfile)

if args.extented_output:
   # Writing the txt file
   with open(out_path+"/results.txt", 'w') as f, \
        open(out_path+"/results_tab.txt", 'w') as tab_file, \
        open(out_path+"/results_table.txt", 'w') as table_file, \
        open(out_path+"/Plasmid_seq.fsa", 'w') as ref_file, \
        open(out_path+"/Hit_in_genome_seq.fsa", 'w') as hit_file:
      tab_file.write("".join(tab_file_lst))
      table_file.write("".join(table_file_lst))
      ref_file.write("".join(ref_file_lst))
      hit_file.write("".join(hit_file_lst))

      # Writing table
      f.write(text_table(headers, rows))
   
      # Writing extended output
      if dbs_with_results:
         f.write("\n\nExtended Output:\n")
         f.write("%s\n"%('-'*80))
         for db_name in dbs_with_results:
            # Txt file alignments
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

