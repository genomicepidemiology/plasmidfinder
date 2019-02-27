#!/usr/bin/env python3
from __future__ import division
from argparse import ArgumentParser
from tabulate import tabulate
from cgecore.blaster import Blaster
from cgecore.cgefinder import CGEFinder
from distutils.spawn import find_executable
import sys, os, time, re, subprocess
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
   rows = [[x if x else empty_replace for x in row] for row in rows]
   # Create table
   table = tabulate(rows, headers, tablefmt='simple').split('\n')
   # Prepare title injection
   width = len(table[0])
   # Switch horisontal line
   table[1] = '*'*(width+2)
   # Update table with title
   table = ("%s\n"*3)%('*'*(width+2), '\n'.join(table), '='*(width+2))
   return table

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

def make_aln(file_handle, json_data, query_aligns, homol_aligns, sbjct_aligns):
   for dbs_info in json_data.values():
      for db_name, db_info in dbs_info.items():
         if isinstance(db_info, str):
            continue
         for gene_id, gene_info in sorted(db_info.items(), key=lambda  x: (x[1]['plasmid'], x[1]['accession'])):  
            seq_name = gene_info["plasmid"] + "_" + gene_info["accession"]
            hit_name = gene_info["hit_id"]

            seqs = ["","",""]
            seqs[0] = sbjct_aligns[db_name][hit_name]    
            seqs[1] = homol_aligns[db_name][hit_name]    
            seqs[2] = query_aligns[db_name][hit_name]

            write_align(seqs, seq_name, file_handle)

def write_align(seq, seq_name, file_handle):
   file_handle.write("# {}".format(seq_name) + "\n")
   sbjct_seq = seq[0]
   homol_seq = seq[1]
   query_seq = seq[2]
   for i in range(0,len(sbjct_seq),60):
      file_handle.write("%-10s\t%s\n"%("template:", sbjct_seq[i:i+60]))
      file_handle.write("%-10s\t%s\n"%("", homol_seq[i:i+60]))
      file_handle.write("%-10s\t%s\n\n"%("query:", query_seq[i:i+60]))


##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()
parser.add_argument("-i", "--infile", dest="infile", help="FASTA or FASTQ input files.", nargs = "+")
parser.add_argument("-o", "--outputPath", dest="outdir",help="Path to blast output", default='.')
parser.add_argument("-tmp", "--tmp_dir", help="Temporary directory for storage of the results from the external software.")
parser.add_argument("-mp", "--methodPath", dest="method_path",help="Path to method to use (kma or blastn)")
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='/database')
parser.add_argument("-d", "--databases", dest="databases",help="Databases chosen to search in - if non is specified all is used")
parser.add_argument("-l", "--mincov", dest="min_cov",help="Minimum coverage", default=0.60)
parser.add_argument("-t", "--threshold", dest="threshold",help="Minimum threshold for identity", default=0.90)
parser.add_argument("-x", "--extented_output",
                    help="Give extented output with allignment files, template and query hits in fasta and\
                          a tab seperated file with allele profile results", action="store_true")
parser.add_argument("-q", "--quiet", action="store_true")

args = parser.parse_args()

##########################################################################
# MAIN
##########################################################################

if args.quiet:
   f = open('/dev/null', 'w')
   sys.stdout = f

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
   db_path = args.db_path

# Check if valid input files are provided
if args.infile is None:
   sys.exit("Input Error: No input file was provided!\n")
elif not os.path.exists(args.infile[0]):
   sys.exit("Input Error: Input file does not exist!\n")
elif len(args.infile) > 1:
   if not os.path.exists(args.infile[1]):
      sys.exit("Input Error: Input file does not exist!\n")
   infile = args.infile
else:   
   infile = args.infile

# Check if valid output directory is provided
if not os.path.exists(args.outdir):
   sys.exit("Input Error: Output dirctory does not exist!\n")
outdir = os.path.abspath(args.outdir)

# Check if valid tmp directory is provided
if args.tmp_dir:
   if not os.path.exists(args.tmp_dir):
      sys.exit("Input Error: Tmp dirctory, {}, does not exist!\n".format(args.tmp_dir))
   else:
      tmp_dir = os.path.abspath(args.tmp_dir)
else:
   tmp_dir = outdir

# Check if databases and config file are correct/correponds
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

      # Check if all db files are present
      for ext in extensions:
         db = "%s/%s.%s"%(db_path, db_prefix, ext)
         if not os.path.exists(db):
            sys.exit(("Input Error: The database file (%s) "
                                 "could not be found!")%(db))
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

species = [",".join(dbs[db]) for db in databases]


# Check file format (fasta, fastq or other format)
file_format = get_file_format(infile)

# Call appropriate method (kma or blastn) based on file format 
if file_format == "fastq":
   if not method_path:
      method_path = "kma"
   if find_executable(method_path) == None:
      sys.exit("No valid path to a kma program was provided. Use the -mp flag to provide the path.")

   # Check the number of files
   if len(infile) == 1:
      infile_1 = infile[0]
      infile_2 = None
   elif len(infile) == 2:
      infile_1 = infile[0]
      infile_2 = infile[1]
   else:
      sys.exit("Only 2 input file accepted for raw read data,\
                if data from more runs is avaliable for the same\
                sample, please concatinate the reads into two files")
    
   sample_name = os.path.basename(sorted(args.infile)[0])
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
   assert len(infile) == 1, "Only one input file accepted for assembled data"
   infile = infile[0]
   method = "blast"

   # Call BLASTn
   method_obj = Blaster(infile, databases, db_path, tmp_dir, 
                        min_cov, threshold, method_path, cut_off=False)
else:
   sys.exit("Input file must be fastq or fasta format, not "+ file_format)

results      = method_obj.results
query_aligns = method_obj.gene_align_query
homo_aligns  = method_obj.gene_align_homo
sbjct_aligns = method_obj.gene_align_sbjct

json_results = dict()

hits = []

for db in results:
   contig_res = {}
   if db == 'excluded':
      continue
   db_name = str(dbs[db][0])
   if db_name not in json_results:
      json_results[db_name] = {}

   if db not in json_results[db_name]: 
      json_results[db_name][db] = {}

   if results[db] == "No hit found":
      json_results[db_name][db] = "No hit found"
   else:
      for contig_id, hit in results[db].items():

         identity = float(hit["perc_ident"])
         coverage = float(hit["perc_coverage"])

         # Skip hits below coverage
         if coverage < (min_cov*100) or identity < (threshold*100):
            continue

         bit_score = identity * coverage

         if contig_id not in contig_res:
            contig_res[contig_id] = []
         contig_res[contig_id].append([hit["query_start"], hit["query_end"], bit_score, hit])

   # Check for overlapping hits, only report the best   
   for contig_id, hit_lsts in contig_res.items():

      hit_lsts.sort(key=lambda x: x[0])
      hits = [hit[3] for hit in hit_lsts]

      # Get information from the fisrt hit found
      current_end = hit_lsts[0][1]
      current_bit_score = hit_lsts[0][2]
    
      # Check if more then one hit was found within the same gene
      for i in range(len(hit_lsts)-1):

         # Save information from next hit
         next_start = hit_lsts[i+1][0]
         next_end = hit_lsts[i+1][1]
         next_bit_score = hit_lsts[i+1][2]

         # Check for overlapping sequences
         # <--------------->
         #            <------------>
         if next_start < current_end:
            # Delete the hit with lowest bit score from the hit list
            #<----->
            #   <-------------------->
            if current_bit_score < next_bit_score:
               del hits[i]
               # reset current end and score
               current_end = next_end
               current_bit_score = next_bit_score
            #<-------------------->
            #                  <----->
            else:
               # Delete next hit, keep current end and bit score
                del hits[i+1]

      for hit in hits:
         header = hit["sbjct_header"]

         tmp = header.split("_")
         try:
            gene = tmp[0]
            note = tmp[2]
            acc = "".join(tmp[3:])
         except IndexError:
            gene = tmp
            note = ""
            acc = ""
         identity = hit["perc_ident"]
         coverage = hit["perc_coverage"]
         sbj_length = hit["sbjct_length"]
         HSP = hit["HSP_length"]
         positions_contig = "%s..%s"%(hit["query_start"], hit["query_end"])
         positions_ref = "%s..%s"%(hit["sbjct_start"], hit["sbjct_end"])
         contig_name = hit["contig_name"]

         # Write JSON results dict
         json_results[db_name][db].update({header:{}})
         json_results[db_name][db][header] = {"plasmid":gene,"identity":round(identity, 2),"HSP_length":HSP,
                                              "template_length":sbj_length,"position_in_ref":positions_ref,
                                              "contig_name":contig_name,"positions_in_contig":positions_contig,
                                              "note":note,"accession":acc, 
                                              "coverage":round(coverage, 2), "hit_id":contig_id}

# Get run info for JSON file
service = os.path.basename(__file__).replace(".py", "")
date = time.strftime("%d.%m.%Y")
time = time.strftime("%H:%M:%S")

# Make JSON output file
data = {service:{}}

userinput = {"filename(s)":args.infile, "method":method,"file_format":file_format}
run_info = {"date":date, "time":time}

data[service]["user_input"] = userinput
data[service]["run_info"] = run_info
data[service]["results"] = json_results

pprint.pprint(data)

# Save json output
result_file = "{}/data.json".format(outdir) 
with open(result_file, "w") as outfile:
   json.dump(data, outfile)

# Getting and writing out the results
header = ["Plasmid", "Identity", "Query / Template length", "Contig", "Position in contig", "Note", "Accession number"]

if args.extented_output:
   # Define extented output 
   table_filename  = "{}/results_tab.tsv".format(outdir)
   query_filename  = "{}/Hit_in_genome_seq.fsa".format(outdir)
   sbjct_filename  = "{}/Plasmid_seqs.fsa".format(outdir)
   result_filename = "{}/results.txt".format(outdir)
   table_file  = open(table_filename, "w")
   query_file  = open(query_filename, "w")
   sbjct_file  = open(sbjct_filename, "w")
   result_file = open(result_filename, "w")

   # Make results file
   result_file.write("{} Results\n\nOrganism(s): {}\n\n".format(service, ",".join(set(species))))

   # Write tsv table
   rows = [["Database"] + header]
   for species, dbs_info in json_results.items():
      for db_name, db_hits in dbs_info.items():
         result_file.write("*"*len("\t".join(header)) + "\n")
         result_file.write(",".join(dbs[db_name]) + "\n")
         db_rows = []

         # Check it hits are found
         if isinstance(db_hits, str):
            content = ['']*len(header)
            content[int(len(header)/2)] = db_hits
            result_file.write(text_table(header, [content]) + "\n")
            continue

         for gene_id, gene_info in sorted(db_hits.items(), key=lambda  x: (x[1]['plasmid'], x[1]['accession'])):
            vir_gene = gene_info["plasmid"]
            identity = str(gene_info["identity"])
            coverage = str(gene_info["coverage"])
            template_HSP = str(gene_info["HSP_length"]) + " / " + str(gene_info["template_length"])
            position_in_ref = gene_info["position_in_ref"]
            position_in_contig = gene_info["positions_in_contig"]
            note = gene_info["note"]
            acc = gene_info["accession"]
            contig_name = gene_info["contig_name"]

            # Add rows to result tables
            db_rows.append([vir_gene, identity, template_HSP, contig_name, position_in_contig, note, acc])
            rows.append([db_name, vir_gene, identity, template_HSP, contig_name, position_in_contig, note, acc])

            # Write query fasta output
            hit_name = gene_info["hit_id"]
            query_seq = query_aligns[db_name][hit_name]
            sbjct_seq = sbjct_aligns[db_name][hit_name] 

            if coverage == 100 and identity == 100:
               match = "PERFECT MATCH"
            else:
               match = "WARNING"
            qry_header = ">{}:{} ID:{}% COV:{}% Best_match:{}\n".format(vir_gene, match, identity, 
                                                                    coverage, gene_id)
            query_file.write(qry_header)
            for i in range(0,len(query_seq),60):
               query_file.write(query_seq[i:i+60] + "\n")

            # Write template fasta output
            sbj_header = ">{}\n".format(gene_id)
            sbjct_file.write(sbj_header)
            for i in range(0,len(sbjct_seq),60):
               sbjct_file.write(sbjct_seq[i:i+60] + "\n")

         # Write db results tables in results file and table file
         result_file.write(text_table(header, db_rows) + "\n")

      result_file.write("\n")

   for row in rows:
      table_file.write("\t".join(row) + "\n")

   # Write allignment output
   result_file.write("\n\nExtended Output:\n\n")
   make_aln(result_file, json_results, query_aligns, homo_aligns, sbjct_aligns)

    # Close all files
   query_file.close()
   sbjct_file.close()
   table_file.close()
   result_file.close()

if args.quiet:
   f.close()
