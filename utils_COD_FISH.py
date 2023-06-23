#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 09:20:13 2021

@author: Ryan Dikdan, Devisi Goel

This is the utility script which has functions for COD_FISH.py which
takes gene names and generates smFISH probes for the gene.
The motivation is to find optimal oligos for targeting RNAs given full transcriptome data.
"""

##############################################################################
#   Imports                                                                  #
##############################################################################
import subprocess
import time
import primer3
import requests
import math
import os, sys
import platform
import multiprocessing as mp
import gzip
import shutil
import psutil
import ensembl_rest
from tqdm import tqdm
import urllib.request   # for ftp downloads

##############################################################################
#  DNA helper functions                                                      #
##############################################################################
def complement(seq):
    complement_dict = {'A':'T','T':'A','C':'G','G':'C'}
    complement_list = []
    for nt in seq.replace('-','').upper():  # Removes -'s from BLAST. This should bring back the actual sequence, ignoring gaps.
        complement_list.append(complement_dict[nt])
    complement_seq = ''.join(complement_list)
    return(complement_seq)

def reverse_complement(seq):
    return(complement(seq)[::-1])


##############################################################################
#  File management functions                                                 #
##############################################################################

def download_file(url_path, output_filename):
    with requests.get(url_path, stream=True) as file_request:
        total_length = int(file_request.headers.get("Content-Length"))
        with tqdm.wrapattr(file_request.raw, "read", total=total_length, desc="") as raw_data:
            with open(output_filename, 'wb') as output_file_handle:
                shutil.copyfileobj(raw_data, output_file_handle)

def gunzip(gz_file):
    # gz Unzip
    with gzip.open(gz_file, 'rb') as file_in:
        # '.'.join(gz_file.split('.')[:-1]) removes the '.gz'
        with open('.'.join(gz_file.split('.')[:-1]), 'wb') as file_out:
            shutil.copyfileobj(file_in, file_out)

##############################################################################
#   Setup of config, species, and input genes                                #
##############################################################################

def make_config_file(num_probes = 50, probe_len = 20, probe_dist = 2, min_tm = 40, hairpin_max_tm = 57, offtarget_calc_method = "Tm", probe_selection_method = "Greedy", filter_repeats = True):
    # Get current directory
    # This loop checks if it's run through pyinstaller and picks the appropriate path
    if getattr(sys, 'frozen', False):
        # If the application is run as a bundle, the PyInstaller bootloader
        # extends the sys module by a flag frozen=True and sets the app 
        # path into variable _MEIPASS'.
        application_path = sys._MEIPASS
    else:
        application_path = os.path.dirname(os.path.abspath(__file__))
    # Changes workign directory to the newly made species directory
    os.chdir(application_path)
    config_file = open("./config.py", "w")
    config_file.write("# This is the configuration file which defines the various parameters for making probes\n\n")

    config_file.write("time_made = \'" + time.asctime(time.localtime(time.time()))+'\'\n')
    config_file.write("num_probes = "+str(num_probes)+" #Number of probes to be created\n")
    config_file.write("probe_len = "+str(probe_len)+" #Length of each probe\n")
    config_file.write("probe_dist = "+str(probe_dist)+" #Distance between each probe\n")
    config_file.write("min_tm = "+str(273+float(min_tm))+" #Minimum melting temperature of probes in Kelvin\n")
    config_file.write("hairpin_max_tm = "+str(273+float(hairpin_max_tm))+" #Maximum melting temperature of the hairpin of probes in Kelvin\n\n")

    config_file.write("offtarget_calc_method = \'" + offtarget_calc_method+'\'\n')
    config_file.write("probe_selection_method = \'" + probe_selection_method +'\'\n')
    
    config_file.write("#Filter repeating nucleotides, AAAAA, GGGG, CCCC, TTTTT\n")
    config_file.write("filter_repeats = "+str(filter_repeats)+"\n")

    config_file.close()

def make_and_check_species_folder(species):
    # Get current directory
    # This loop checks if it's run through pyinstaller and picks the appropriate path
    if getattr(sys, 'frozen', False):
        # If the application is run as a bundle, the PyInstaller bootloader
        # extends the sys module by a flag frozen=True and sets the app 
        # path into variable _MEIPASS'.
        application_path = sys._MEIPASS
    else:
        application_path = os.path.dirname(os.path.abspath(__file__))
    # Changes workign directory to the correct one
    os.chdir(application_path)

    # Making the directory for the transcriptome files
    species_dir = "species." + species + "/"
    if not os.path.exists(application_path+'/'+species_dir):
        print('Making the species specific directory for', species, '\n')
        os.mkdir(application_path+'/'+species_dir)
    else:
        print('Species specific directory for', species, 'already made\n')
    sys.path.append(application_path+'/'+species_dir)

    # Changes workign directory to the newly made species directory
    os.chdir(species_dir)

    
    # Download appropriate files through ensembl using curl. Not at all intuitive.
    # Uses curl to pull all files in the species cDNA folder.
    # Then find the filename that has the appropriate name (containing cdna.fa.gz in this case)
    # then downloading it.

    raw_cDNA_site_data = requests.get('http://ftp.ensembl.org/pub/current_fasta/'+species+'/cdna/')
    cDNA_site_data = raw_cDNA_site_data.text.split("\"")
    for entry in cDNA_site_data:
        if 'cdna.all.fa.gz' in entry and '</a>' not in entry:
            cDNA_file_name = entry
    if len(cDNA_file_name) < 3:
        raise(Exception(
            'Couldn\'t find the appropriate files on ensemble. Please contact support'))
    if not os.path.exists('cDNA.fa'):
        print("Downloading the cDNA file\n")
        download_file('http://ftp.ensembl.org/pub/current_fasta/'+species+'/cdna/'+cDNA_file_name, 'cDNA.fa.gz')
        
        # with requests.get('http://ftp.ensembl.org/pub/current_fasta/'+species+'/cdna/'+cDNA_file_name, stream=True) as cDNA_request:
        #     with open('cDNA.fa.gz', 'wb') as cDNA_gz_file:
        #         shutil.copyfileobj(cDNA_request.raw, cDNA_gz_file)
        
        gunzip('cDNA.fa.gz')

    else:
        print("cDNA file already downloaded")

    raw_ncRNA_site_data = requests.get('http://ftp.ensembl.org/pub/current_fasta/'+species+'/ncrna/')
    ncRNA_site_data = raw_ncRNA_site_data.text.split("\"")
    for entry in ncRNA_site_data:
        if 'ncrna.fa.gz' in entry and '</a>' not in entry:
            ncRNA_file_name = entry
    if len(ncRNA_file_name) < 3:
        raise(Exception(
            'Couldn\'t find the appropriate files on ensemble. Please contact support'))
    if not os.path.exists('ncRNA.fa'):
        print("Downloading the ncRNA file\n")
        
        download_file('http://ftp.ensembl.org/pub/current_fasta/'+species+'/ncrna/'+ncRNA_file_name, 'ncRNA.fa.gz')
        
        # with requests.get('http://ftp.ensembl.org/pub/current_fasta/'+species+'/ncrna/'+ncRNA_file_name, stream=True) as ncRNA_request:
        #     with open('ncRNA.fa.gz', 'wb') as ncRNA_gz_file:
        #         shutil.copyfileobj(ncRNA_request.raw, ncRNA_gz_file)

        gunzip('ncRNA.fa.gz')
        
    else:
        print("ncRNA file already downloaded")

    # Combining the cDNA and ncRNA files for simplicity and then extracting it
    if not os.path.exists('transcriptome.fa'):
        print('Combining cDNA and ncRNA files file')
        
        # Using a weird shutil thing to do it very quickly
        with open('transcriptome.fa','wb') as transcriptome_file:
            for fasta_file in ['cDNA.fa','ncRNA.fa']:
                with open(fasta_file,'rb') as current_fasta:
                    shutil.copyfileobj(current_fasta, transcriptome_file)
    else:
        print('\nUsing already made transcriptome file.\n')

    # Using the generated transcriptome file to make the blast database for alignment of probe sequences
    if not os.path.exists('blast_db.nhr'):
        make_blast_db('transcriptome.fa')
    else:
        print('\nUsing already made blast database files.\n')

    print("The blast database has successfully been made")

    print("The blast database is stored at blast_db")
    print("The transcriptome file is stored at transcriptome.fa")

    os.chdir('../')

def input_to_transcript_dictionary(transcriptome_dict, species,
                                   transcripts_file = None, gene_symbol_file = None, gene_symbols = None, ensembl_gene_ids = None, ensembl_transcript_ids = None,
                                   transcript_targets = {}, ensembl_gene_id_list = [], gene_symbol_list = [], ensembl_transcript_id_list = []):
        # Add user inputed sequences
    if (transcripts_file):
        transcript_targets = load_target_transcripts(transcripts_file)

    # Add gene symbols from the command line
    if (gene_symbols):
        gene_symbol_list += gene_symbols

    # Add gene symbols from a file
    if (gene_symbol_file):
        for line in open(gene_symbol_file):
            gene_symbol = line.strip()
            gene_symbol_list.append(gene_symbol)

    # Add gene IDs from the command line
    if (ensembl_gene_ids):
        for ensembl_gene_id in list(ensembl_gene_ids):
            ensembl_gene_id_list.append(ensembl_gene_id)

    # If gene symbols or gene IDs were input then use pybiomart to find the primary transcript IDs
    if len(gene_symbol_list+ensembl_gene_id_list) > 0:
        primary_ensembl_transcript_ids = load_gene_information(
            gene_symbol_list, ensembl_gene_id_list, species)
        for ensembl_transcript_id in primary_ensembl_transcript_ids:
            ensembl_transcript_id_list.append(ensembl_transcript_id)

    # If transcript IDs were input then they are simply added to the transcript list
    if (ensembl_transcript_ids):
        for ensembl_transcript_id in list(ensembl_transcript_ids):
            ensembl_transcript_id_list.append(ensembl_transcript_id)

    # If there are transcript IDs in the transcript list then it loads them from the dictionary
    if len(ensembl_transcript_id_list) > 0:
        for transcript_id in ensembl_transcript_id_list:
            transcript_targets[transcript_id] = transcriptome_dict[transcript_id]

    # If there were no inputs error out here
    elif len(transcript_targets) == 0:
        raise Exception('No genes input, please make sure to use the -g, -fa, -gid, -tid, or -gf \
                        (--gene_symbols, --fasta_file, --ensembl_gene_ids, --ensembl_transcript_ids, or --gene_symbol_file) options')
    return(transcript_targets)


##############################################################################
#   Loading fastas and finding primary transcripts                           #
##############################################################################

def load_transcriptome_dict(transcriptome_file):

    print('\nLoading the transcriptome into a python dictionary')
    
    transcriptome_dict = {}
    gene_info_array = []

    for transcriptome_line in open(transcriptome_file):
        transcriptome_line = transcriptome_line.strip()
        if transcriptome_line[0] == ">":
            transcript_full_name = transcriptome_line[1:]
            transcript_id = transcript_full_name.split(' ')[0].split(".")[0]
            
            gene_symbol = transcript_full_name.partition('gene_symbol:')[2].split(' ')[0]
            if gene_symbol == '':
                gene_symbol = 'not_listed'
            transcript_id = transcript_full_name.split(' ')[0].split(".")[0]
            gene_id = transcript_full_name.partition('gene:')[2].split('.')[0]
            gene_biotype = transcript_full_name.partition('gene_biotype:')[2].split(' ')[0]

            gene_info_array.append([gene_symbol, transcript_id, gene_id, gene_biotype])
            
            transcriptome_dict[transcript_id] = ""      # Makes an entry using the name, where the value is blank
        else:
                                                                # Adds each line to the previously identified ID's sequence
            transcriptome_dict[transcript_id] = transcriptome_dict[transcript_id] + transcriptome_line
    
    # Removing pandas since it makes the compiled version huge even though I'm only using a little bit of it.
    
    # This filters it to make a list of Transcript IDs for rRNAs
    rRNA_list = []
    for gene_row in gene_info_array:
        if gene_row[3] == 'rRNA':
            rRNA_list.append(gene_row[1])
    #
    return(transcriptome_dict, gene_info_array, rRNA_list)

# Loads a user inputted fasta file to generate probes against it
def load_target_transcripts(input_fasta):
    # Parsing, converting fasta to dictionary
    # Fasta should have sequence names of the format:
    # >transcript_id for example:
    #
    # >ENST00000681605
    # AATTGTCATACGACTTGCAGTGAGCGTCAGGAGCACGTCCAGGAACTCCTCAGCAGCGCC
    # TCCTTCAGCTCCACAGCCAGACGCCCTCAGACAGCAAAGCCTACCCCCGCGCCGCGCCCT
    # ...
    # This is especially useful if you want to target a specific 'non-canonical' transcript
    # This info is used so that we can exclude binding penalties to the target gene
    #
    # If you would like to ignore this you can just list an ID like this >XXXXXXXXXXXXXXXXXXX
    
    target_dict = {}
    for line in open(input_fasta):
        line = line.strip()
        if line[0] == ">":
            transcript_name = line[1:]
            target_dict[transcript_name] = ''
        else:
            target_dict[transcript_name] = target_dict[transcript_name] + line.upper()
    return(target_dict)

# From ensembl REST API via ensembl_rest
def load_gene_information(gene_symbols_list, gene_id_list, species):
    
    print('\nLoading transcript IDs for genes by symbol and ID through ensembl_rest')
    
    primary_transcript_list = []
    
    # Loads canonical transcripts from gene symbols, ie PTGS2, NOS2
    for gene_symbol in gene_symbols_list:
        primary_transcript_list.append(ensembl_rest.symbol_lookup(species=species,symbol=gene_symbol)['canonical_transcript'].split('.')[0])
    
    # Loads canonical transctips from gene IDs, ie ENSG00000073756, ENSG00000007171
    for gene_ID in gene_id_list:
        primary_transcript_list.append(ensembl_rest.lookup(species=species,symbol=gene_ID)['canonical_transcript'].split('.')[0])

    return(primary_transcript_list)


##############################################################################
# Create Probe Possibilities, filter by Tm and Tm_hairpin                    #
##############################################################################

# Takes the config parameters and filters out probes
def create_and_filter_probe_possibilites(target_transcript_seq, probe_len, Tm, Tm_hairpin, num_probes, filter_repeats = True):
    
    print("Making and filtering probe possibilities by Tm, repeats, and hairpins\n")

    probe_seq_dict = {}
    filtered_by_tm = 0
    filtered_by_repeat = 0
    filtered_by_hairpin = 0

    # Make all probe possibilities
    candidate_probe_list = []
    for location in range(0,len(target_transcript_seq)-probe_len + 1):
        probe_region = target_transcript_seq[location:location+probe_len]
        probe_seq = reverse_complement(probe_region)
        probe_name = "probe_" + str(location+1)
        candidate_probe_list.append((probe_name, probe_seq))

    # Filter out if Tm is too low
    Tm_filtered_probes_list = []
    for candidate_probe in candidate_probe_list:
        probe_name = candidate_probe[0]
        probe_seq = candidate_probe[1]
        probe_tm = 273 + primer3.bindings.calc_tm(probe_seq)
        if float(probe_tm) > Tm:
            Tm_filtered_probes_list.append((probe_name, probe_seq, probe_tm))
        else:
            filtered_by_tm += 1
    if len(Tm_filtered_probes_list) < num_probes:
        raise Exception("Error: Too many probes filtered out due to low probe binding Tm, \
                        consider lowering the minimum Tm in the config.py file and rerunning.")

    # Filter out if repetitive sequences
    Tm_rep_filtered_probes = []
    if filter_repeats == True:
        for candidate_probe in Tm_filtered_probes_list:
            probe_name = candidate_probe[0]
            probe_seq = candidate_probe[1] 
            probe_tm = candidate_probe[2]
            if ('AAAAA' or 'TTTTT' or 'GGGG' or 'CCCC') not in probe_seq:
                Tm_rep_filtered_probes.append((probe_name,probe_seq,probe_tm))
            else:
                filtered_by_repeat += 1
    else:
        Tm_rep_filtered_probes = Tm_filtered_probes_list
    if len(Tm_rep_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to repetitive nucleotides (ie. AAAAA or GGGG), \
                        consider changing \"filter_repeats = True\" to \"filter_repeats = False\" in the config.py file and rerunning.")

    # Filter out if high hairpin Tm
    Tm_rep_hairpin_filtered_probes = []
    for candidate_probe in Tm_rep_filtered_probes:
        probe_name = candidate_probe[0]
        probe_seq = candidate_probe[1]
        probe_tm = candidate_probe[2]
        probe_hairpin_tm = primer3.bindings.calc_hairpin(probe_seq).tm + 273
        if float(probe_hairpin_tm) < Tm_hairpin:
            Tm_rep_hairpin_filtered_probes.append((probe_name,probe_seq,probe_tm,probe_hairpin_tm))
            probe_seq_dict[probe_name] = probe_seq
        else:
            filtered_by_hairpin += 1
    if len(Tm_rep_hairpin_filtered_probes) < num_probes:
        raise Exception("Error: Too many probes filtered out due to repetitive nucleotides (ie. AAAAA or GGGG), \
                        consider changing \"filter_repeats = True\" to \"filter_repeats = False\" in the config.py file and rerunning.")
    
    print("\nTm filtered: "+str(filtered_by_tm))
    print("Repeat filtered: "+str(filtered_by_repeat))
    print("Hairpin filtered: "+str(filtered_by_hairpin)+'\n')
        
    return(Tm_rep_hairpin_filtered_probes, probe_seq_dict)


##############################################################################
#  Scoring functions                                                         #
##############################################################################

# This function is used by .sort to sort the tuples based on their score
# All it does is take a tuple and return the 2nd value which is then used to sort
def sort_by_score(tuple):

    return(tuple[1])

def compute_offtarget_scores(alignment_data, target_ensembl_transcript_id, rRNA_list):
    print('Calculating probe candidate off-target scores using blast scores\n')
    offtarget_match_scores = {}

    for aln in alignment_data:
        probe_name = aln[0]
        # Pull ensemble transcript ID from name
        aligned_ensembl_transcript_id = aln[1].split(".")[0]
        #if aligned_ensembl_transcript_id == "*":
            #continue
            # This skips if the probe didn't align to anything
            # I don't believe that this happens with blast
        if aligned_ensembl_transcript_id == target_ensembl_transcript_id:
            continue
            # Doesn't add anything to the probe score if the alignment is an on target alignment

        match_score = int(aln[4])        # This is the score from blast which is the Raw Score

        # Significantly punishes for rRNA alignment. This should make it so these probes aren't selected at all
        if aligned_ensembl_transcript_id in rRNA_list:
            match_score = match_score + 10000

        if probe_name not in offtarget_match_scores:
            offtarget_match_scores[probe_name] = 0
            # Adds entry to offtarget scores if probe isn't already in it
        
        # Add the score to the off target score for that probe
        offtarget_match_scores[probe_name] = offtarget_match_scores[probe_name] + match_score

    # Converts dictionary to an array of tuples then sorts it based on the scores
    offtarget_match_scores_tuple_list = list(offtarget_match_scores.items())
    offtarget_match_scores_tuple_list.sort(key=sort_by_score)

    return(offtarget_match_scores_tuple_list)


def compute_offtarget_tm_for_mp(aln):

    probe_name = aln[0]
    aligned_ensembl_transcript_id = aln[1].split(".")[0]
    probe_seq = aln[2]
    aligned_transcript_seq_fragment = reverse_complement(aln[3])

    aln_tm = (primer3.bindings.calc_heterodimer(probe_seq, aligned_transcript_seq_fragment).tm + 273)/(primer3.bindings.calc_tm(probe_seq) + 273) + 1
    # Punished 1 per alignment as well as 0-1 based on the alignment 
    # (so for a perfect off target alignment penalty of 2 and for a very weak off target alignment penalty of 1)
    
    return(probe_name, aln_tm, aligned_ensembl_transcript_id)

def multithread_offtarget_tm(alignment_data, target_ensembl_id, rRNA_list):
    
    print('Calculating probe candidate off-target scores using calculated alignment Tm\'s using multithreading\n')
    # Function which just calculates Tm penalty which is used by the multithreading Pool

    
    output_probe_score = {}
    # Here's the multithreading which takes the compute_offtarget_tm_for_mp function above and applies it to every alignment
    with mp.Pool() as executor:
        tm_scores = executor.map(compute_offtarget_tm_for_mp, alignment_data)
    # Then we take the resulting tm_scores variable and combine the alignments per probe 
    # and using some conditional statements on them
    #
    # Tried to make faster by filtering out on targets, but then it requires 
    # 2 loops through and ends up being a tad slower.
    for score_result_unit in tm_scores:
        probe_name = score_result_unit[0]
        probe_score = score_result_unit[1]
        aligned_ensembl_id = score_result_unit[2]
        
        # Doesn't add anything to the probe score if the alignment is an on target alignment
        if aligned_ensembl_id == target_ensembl_id:
            continue
        
        # Significantly punishes if aligned to an rRNA
        if aligned_ensembl_id in rRNA_list:
            probe_score = probe_score + (probe_score + 10000)
        
        # Adds probe to output if not there already
        if probe_name not in output_probe_score.keys():
            output_probe_score[probe_name] = probe_score
        
        # Adds score if already there
        elif output_probe_score[probe_name]:
            output_probe_score[probe_name] = output_probe_score[probe_name] + probe_score

    offtarget_tuple_list = []
    # Add the values to a list so we can sort it (dictionaries are unsorted)
    for entry in output_probe_score:
        offtarget_tuple_list.append((entry,output_probe_score[entry]))
    offtarget_tuple_list.sort(key=sort_by_score)
    
    return(offtarget_tuple_list)

##############################################################################
#  Blast functions                                                         #
##############################################################################

def check_blast_download_install():
    # This checks if blastn in a $PATH variable, which pretty much checks if it's installed and usable.
    if not shutil.which('blastn'):
        print('Blast isn\'t installed, to download and install it. Please follow the installation instructions.')
        print('This program or the computer might need to be restarted for blast to function properly after installation.')
        
        # If run on windows, will download and install the windows version of blast
        if os.name == 'nt':
            urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-win64.exe', 'blast-2.13.0+-win64_installer.exe')
            blast_install = subprocess.call(['blast-2.13.0+-win64_installer.exe'], shell=True)  # Apparently the shell=True part is important for some reason
                                                                                                # that I think has to do with permissions to run .exe files.
            # If the installer doesn't finish it'll throw an error which will be picked up in the blast_install variable. It is 0 if successful.
            if blast_install == 1:
                raise(Exception('The blast installation didn\'t appear to work. Please rerun the script to install blast, or manually download it. \
                                Here\'s the download website: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'))
            
        # Mac version
        if platform.system() == 'Darwin':
            # Should change this to be more robust in case the blast dmg changes.
            urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+.dmg', 'ncbi-blast-2.13.0+.dmg')
            # Mount the dmg
            subprocess.call(['hdiutil', 'attach', '-mountpoint', '/Volumes/blast_installer', 'ncbi-blast-2.13.0+.dmg'])

            # Copy the uninstaller and run the installer
            #subprocess.call(['cp','/Volumes/blast_installer/uninstall_ncbi_blast.zip','blast_uninstaller.zip'])
            # Uninstaller doesn't work on more recent versions of Mac I think.
            print('Please follow the BLAST installer instructions. Installation of BLAST is required for COD-FISH to function.')
            subprocess.call(['open','/Volumes/blast_installer/ncbi-blast-2.13.0+.pkg'])
            
            # To wait for it to complete we will wait based on its name with this function
            def wait_for_process_by_name(process_name):
                """
                Wait for a process with the given name to finish.
                """
                while True:
                    for proc in psutil.process_iter(['pid', 'name']):
                        if proc.name() == process_name:
                            # The process is still running, sleep for a short time and check again
                            time.sleep(0.1)
                            break
                    else:
                        # The process has finished, break out of the loop and continue with the rest of the script
                        break
                        
            wait_for_process_by_name("Installer")
            
            # If the user didn't select removing the volume and dmg, these lines will do that. The errors they
            # throw if the user did select to remove them won't cause any issues
            subprocess.call(['hdiutil', 'detach', '/Volumes/blast_installer'])
            # Remove the dmg file
            subprocess.call(['rm','ncbi-blast-2.13.0+.dmg'])
                
    

def make_blast_db(fasta_file):
    check_blast_download_install()
    print("\nBlastn will now generate the database with which to align the candidate probe sequences.")
    
    # If running windows will tell it not to open another terminal, otherwise don't
    # By default other OS don't open new terminals
    if os.name == 'nt':
        blast_db_run = subprocess.run(['makeblastdb',
                                '-in', fasta_file,
                                '-dbtype', 'nucl',
                                '-out','blast_db'],
                                capture_output=True, 
                                creationflags = subprocess.CREATE_NO_WINDOW)
    # Mac
    if platform.system() == 'Darwin':
        blast_db_run = subprocess.run(['makeblastdb',
                                '-in', fasta_file,
                                '-dbtype', 'nucl',
                                '-out','blast_db'],
                                capture_output=True)
    # Build the blast database through subprocess.run
    
    # Return exception if didn't complete normally
    if blast_db_run.returncode != 0:
        raise Exception('Building the database didn\'t work as planned. Make sure the makeblastdb version you have is up to date and all the necessary files are present and in the right folders.')
    return(blast_db_run)  

def align_blast(probe_fasta_filename, transcriptome_db, probe_len):
    check_blast_download_install()
    print('Aligning probe candidates to transcriptome using blast\n')
    # To understand how the scoring and alignment works look at this site: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/
    # In summary, for the raw score 'score' in task blastn-short (which we're doing) the reward is 1 for each match, and -3 for each mismatch. 
    # If a gap is present then it penalizes -5 for gap opening and -2 for extending the gap.
    #
    # Here subprocess runs the actual alignment on the minus strand, with threads equal to what the compute has outputting into a useful format for this program.
    
    # If length of probes is too high will change from blastn-short to blastn. Cutoff is 30nt, as specified on BLAST website: https://www.ncbi.nlm.nih.gov/books/NBK279684/
    if probe_len >= 30:
        task = 'blastn'
    else:
        task = 'blastn-short'
        
    # If running windows will tell it not to open another terminal, otherwise don't
    # By default other OS don't open new terminals
    
    if os.name == 'nt':
        blast_run = subprocess.run(['blastn',
                                '-task', task,
                                '-db', transcriptome_db,
                                '-query', probe_fasta_filename, 
                                '-num_threads', str(os.cpu_count()), 
                                '-strand', 'minus',
                                
                                '-outfmt','6 qacc sacc qseq sseq score'],
                                capture_output=True, 
                                creationflags = subprocess.CREATE_NO_WINDOW)
    if platform.system() == 'Darwin':
        blast_run = subprocess.run(['blastn',
                                '-task', task,
                                '-db', transcriptome_db,
                                '-query', probe_fasta_filename, 
                                '-num_threads', str(os.cpu_count()), 
                                '-strand', 'minus',
                                
                                '-outfmt','6 qacc sacc qseq sseq score'],
                                capture_output=True)

    raw_data = blast_run.stdout.decode("utf-8").split('\n')    # Converts to string then to a list where each element is a line
    raw_data.remove('')
    alignment_data = []
    for line in raw_data:
        alignment_data.append(line.split('\t'))               # Splits each alignment line at the tabs then adds it to the alignment_data list
    return(alignment_data, blast_run)


##############################################################################
# Set Selection Method                                                       #
##############################################################################

# Greedy
def select_nonoverlapping_probes(ordered_offtarget_scores_tuple_list, probe_seq_dict, probe_len, probe_dist, probe_num, min_tm):
    # Takes sorted probe tuple list
    nucleotide_idx_covered_already = []
    final_probe_set = []
    total_score = 0.0
    
    for probe_tuple in ordered_offtarget_scores_tuple_list:
        probe_name = probe_tuple[0]
        probe_idx = int(probe_name.split('_')[1])
        probe_score = probe_tuple[1]
        probe_seq = probe_seq_dict[probe_name]

        if probe_idx in nucleotide_idx_covered_already or probe_idx+probe_len in nucleotide_idx_covered_already:
            # Won't add probe if already covered by a probe
            continue
        
        # Test for probe heterodimer binding
        heterodimer_tms = []
        for probe_tuple2 in final_probe_set:
            heterodimer_tms.append(primer3.bindings.calc_heterodimer(probe_seq,probe_tuple2[2]).tm)
        
        for tm in heterodimer_tms:
            if tm > min_tm:
                print('skipping once here')
                continue
                
        final_probe_set.append((probe_name, probe_score, probe_seq))
        total_score += probe_score
        nucleotide_idx_covered_already.extend(range(probe_idx-probe_dist,probe_idx+probe_len+probe_dist))

        if len(final_probe_set) == probe_num:
            final_probe_set.sort(key=sort_by_score)
            return(final_probe_set)
            
    # If couldn't find enough probes
    print("Couldn't find the specified number of probes:", str(probe_num), ". Could only find, ", str(len(final_probe_set))," The program will output the number of probes that it can generate. \
          The transcript may be too small, or the advanced filtering conditions may be too strict. You might be able to decreasing how strict filtering is to get more.")
    return(final_probe_set)


##############################################################################
#  Save probes                                                               #
##############################################################################

def save_probe_file(filtered_probe_array, temp_probe_fasta_filename):

    # Makes the probe file for alignment with blast
    probe_file = open(temp_probe_fasta_filename, "w")
    # Then for each probe outputs it in fasta format
    for probe_tuple in filtered_probe_array:
        probe_name = probe_tuple[0]
        probe_seq = probe_tuple[1]
        probe_file.write(">" + probe_name + "\n")
        probe_file.write(probe_seq + "\n")
    probe_file.close()
    return(temp_probe_fasta_filename)

def save_scores_file(output_folder,transcript, target_gene_symbol, offtarget_scores_list, probe_seq_dict):
    
    # Saves the offtarget scores as a tsv and saves it
    # Makes filename
    probe_scores_file_name = output_folder + target_gene_symbol + '_' + transcript + "_probe_scores.tsv"
    # Opens file
    probe_scores_file = open(probe_scores_file_name, 'w')
    # Writes the data
    for probe in offtarget_scores_list:
        probe_name = probe[0]
        probe_seq = probe_seq_dict[probe_name]
        probe_score = str(probe[1])
        probe_scores_file.write(probe_name + "\t" + probe_score + "\t" + probe_seq + "\n")
    probe_scores_file.close()

def save_output_file(probe_list, output_dir, target_gene_symbol, target_gene_id, transcript, probe_seq_dict, format):
    
    # Generate output filename
    probe_filename = output_dir + target_gene_symbol+"_"+target_gene_id+"_"+transcript + "_probes"
    print("Probe selection is complete. Writing probes to "+ probe_filename)
    # !!! Keep track of how many probes have been written. Might fix this later.
    probes_added = 0
    
    # Exports in the format of choice, either through python .write or through pandas.
    if format == 'fasta':
        probe_file = open(probe_filename+'.fa', "w")
        for probe in probe_list:
            probe_name = probe[0]
            probe_file.write(">" + probe_name + "\n")
            probe_seq = probe_seq_dict[probe_name]
            probe_file.write(probe_seq + "\n")
        probe_file.close()
        
    if format == 'csv':
        probe_file = open(probe_filename+'.csv', "w")
        for probe in probe_list:
            probe_name = probe[0]
            probe_file.write( probe_name + ",")
            probe_seq = probe_seq_dict[probe_name]
            probe_file.write(probe_seq + "\n")
        probe_file.close()
        
    if format == 'tsv':
        probe_file = open(probe_filename+'.tsv', "w")
        for probe in probe_list:
            probe_name = probe[0]
            probe_file.write( probe_name + "\t")
            probe_seq = probe_seq_dict[probe_name]
            probe_file.write(probe_seq + "\n")
        probe_file.close()
        
    if format == 'scores_csv':
        probe_file = open(probe_filename+'.csv', "w")
        for probe in probe_list:
            probe_name = probe[0]
            probe_file.write( probe_name + ",")
            probe_seq = probe_seq_dict[probe_name][0]   # This is wrong 
            probe_file.write(probe_seq + ",")
            probe_score = probe_seq_dict[probe_name][1] # This needs fixing
            probe_file.write(str(probe_score) + "\n")
        probe_file.close()
        
    '''if format == 'excel':
        probe_list_df_format = []
        for probe in probe_list:
            probe_name = probe[0]
            probe_seq = filtered_probe_tm_dict[probe_name][0]
            probe_list_df_format.append([probe_name,probe_seq])
            probes_added += 1
            if probes_added >= num_probes: 
                break
        # Workbook
        wb = openpyxl.Workbook()
        # Worksheet
        ws = wb.active
        ws.title = target_gene_symbol
        for probe_entry in probe_list_df_format:
            ws.append(probe_entry)
        
        wb.save(filename = probe_filename+'.xlsx')'''
        
