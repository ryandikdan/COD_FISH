#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 09:20:13 2021

@author: Ryan Dikdan, Devisi Goel

This is a script which takes gene names and generates smFISH probes for the gene.
The motivation is to find optimal oligos for targeting RNAs given full transcriptome data.
"""
import shutil
import sys
import os
import argparse
import ensembl_rest


import utils_COD_FISH as utils

# Previously, this script could only be run on the command line 
# Then this entire script was added to a main function, which will be run if the script if 
# run on the command line like `python COD_FISH.py` but now since it's part of the main
# function, this main function can be called by other scripts, such as gui_COD_FISH.py

# Look at the bottom to see how this script runs by itself

def main(transcript_targets = {}, species = None, transcriptome_dict = {}, gene_info_array = None, rRNA_list = None):
    # The inputs are what will be done by the Graphical User Interface when this script isn't called directly
    
    # Changes working directory to where the script is so everything lines up well and the directory creation doesn't cause problems
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Check if blast is installed and if not, install it
    utils.check_blast_download_install()

    ##############################################################################
    #  Command Line Interface (CLI) Input                                        #
    ##############################################################################
    
    # Makes parser to collect arguements
    parser = argparse.ArgumentParser(
        description='Generating smFISH probes from given mRNA or gene symbols')

    # Defines all arguements to collect

    # Species
    # Makes the argument required if this script is called directly.
    if __name__ == "__main__":
        parser.add_argument('-s', '--species', metavar='Species', type=str,
                        help='write down the species you want to use (required)',required=True)
    # Otherwise, it isn't required since the species will be input via the main functions inputs
    # and not as an argument
    else:
        parser.add_argument('-s', '--species', metavar='Species', type=str,
                        help='write down the species you want to use (required)',required=False)

    # Genes
    parser.add_argument('-fa', '--transcripts_file',
                        type=str, help='write down all target mRNAs in a fasta file, where the name of each sequence is the ensembl gene ID')
    parser.add_argument('-gid', '--ensembl_gene_ids',  type=str,
                        help='list any number of stable ensembl gene IDs (no decimal) separated by a space (case sensitive!)', nargs='+')
    parser.add_argument('-tid', '--ensembl_transcript_ids',  type=str,
                        help='list any number of stable ensembl transcript IDs (no decimal) separated by a space (case sensitive!)', nargs='+')
    parser.add_argument('-g', '--gene_symbols',  type=str,
                        help='list any number of gene names separated by a space (case sensitive!)', nargs='+')
    parser.add_argument('-gf', '--gene_symbol_file',  type=str,
                        help='text file containing genes to process; each gene symbol should be on a separate line')

    # Others
    parser.add_argument('-a', '--advanced', action='store_true',
                        help='If used will let the user set the probe parameters')
    parser.add_argument('-o', '--output_folder', type=str,
                        help='If specified will let user chose where to save probe output files')
    parser.add_argument('-of', '--output_format', type=str,
                        help='If specified will let user chose format file type, options: fasta (default), csv, tsv')


    # Variable which collects the arguements passed in
    args = parser.parse_args()

    ##############################################################################
    #  Setup output format and folder                                            #
    ##############################################################################

    # Setting output format for probes 
    available_formats = ['fasta','tsv','csv']
    # If input, then set it to the input (tries to figure out from other abbreviations too)
    if args.output_format:
        output_format = args.output_format
        if output_format == 'fsa' or output_format == 'fa':
            output_format = 'fasta'
        if output_format == 'xl' or output_format == 'exl' or output_format == 'xls' or output_format == 'xlsx':
            output_format = 'csv'
    else:
        # If not input then set output to fasta
        output_format = 'fasta'

    # If wrong one selected then raise Exception
    if output_format not in available_formats:
        raise Exception('Invalid output format chosen. Please pick from the following formats: '+str(available_formats))

    # Setting output folder for probes
    if args.output_folder:
        output_folder = args.output_folder + '/'
    else:
        output_folder = './output/'
        

    #########################################################################
    # Load species information and downloads and 
    # processes the necessary files
    #########################################################################
    
    # Species selection settings
    
    # Uses ensembl REST API to pull the available species
    ensembl_species = [species_dict['name'] for species_dict in ensembl_rest.species()['species']]
    
    # To allow pulling species by display names (everything converted to lower case)
    brief_to_ensembl_dict = {}
    for species_dict in ensembl_rest.species()['species']:
        brief_to_ensembl_dict[species_dict['display_name'].lower()] = species_dict['name']
    
    # Checks if species was defined in GUI first
    if not species:
        
        # Checks if the input species is in the display names
        # and if so makes the species variable equal to the corresponding
        # latin species name 
        if args.species in brief_to_ensembl_dict:
            species = brief_to_ensembl_dict[args.species]
        # If input name is in the latin species name, sets species to that
        elif args.species in ensembl_species:
            species = args.species

        else:
            raise Exception("Error\nThe species that you have entered,", args.species, "cannot be found. \
                            You can choose any species in the main ensembl website in the following formats:",
                            ensembl_species[0:5],"\nor \n",brief_to_ensembl_dict.keys()[0:5])

        # Big function from utilities which will make the species folder, download the transcriptome files,
        # combine them, and make the blast database
        utils.make_and_check_species_folder(species)


    #########################################################################
    # Checks advanced flag and if present, will remove current config.py and
    # replace it with the one that is made 
    #########################################################################

    if args.advanced:
        if os.path.exists("config.py"):
            os.remove("config.py")
        num_probes = int(input(
            'Please enter the number of probes you would like to generate per transcript (default 50)\n') or '50')
        probe_len = int(
            input('Please enter the length of the probes to be made (default 20)\n') or '20')
        probe_dist = int(
            input('Please enter the distance between probes (default 2)\n') or '2')
        min_tm = int(input(
            'Please enter the minimum melting temperature for each probe in Celsius (default 40 Celsius)\n') or '40')
        hairpin_max_tm = int(input(
            'Please enter the maximum melting temperature for probe hairpins in Celsius (default 57 Celsius)\n') or '57')
        offtarget_calc_method = str(input(
            "Please enter Tm or Alignment_score to indicate which method you would like to use to evaluate your probes?\n") or 'Tm')
        probe_selection_method = str(input(
            "Please enter Dynamic or Greedy to indicate which probe selection selection you would like to use.\n") or 'Greedy')
        filter_repeats = str(input(
            "Please enter True or False for whether or not probe candidates with nucleotide repeats should be filtered out.\n") or 'True')

        utils.make_config_file(num_probes, probe_len, probe_dist, min_tm, hairpin_max_tm,
                            offtarget_calc_method, probe_selection_method, filter_repeats)
            
        # Move the working directory back to the main folder
        sys.path.insert(0, os.getcwd())
        # Import the new config file values now
        import config
        print('\nUsing config file made at: '+config.time_made)
        print('with the following settings: ')

    ##########################################################################
    # If not then make the default config.py file and import it
    ##########################################################################

    else:
        # If config.py doesn't exist, will make it default then import it.
        if not os.path.exists("config.py"):
            # Uses defaults which are defined in utils and listed here:
            # num_probes=50, probe_len=20, probe_dist=2, min_tm=40, hairpin_max_tm=57, offtarget_calc_method='Tm', probe_selection_method="Greedy", filter_repeats=True
            utils.make_config_file()

        # Move the working directory back to the main folder
        sys.path.insert(0, os.getcwd())
        # Import the new config file values now
        import config
        print('\nUsing config file made at: '+config.time_made)
        print('with the following settings: ')

    ##############################################################################
    #  Configuration Parameters                                                  #
    ##############################################################################
    # Load the config variables
    num_probes = config.num_probes
    probe_len = config.probe_len
    probe_dist = config.probe_dist
    min_tm = config.min_tm
    hairpin_max_tm = config.hairpin_max_tm
    
    offtarget_calc_method = config.offtarget_calc_method
    probe_selection_method = config.probe_selection_method
    
    filter_repeats = config.filter_repeats

    # Check that the values passed in are valid
    possible_offtarget_calc_method = ['Alignment_score','Tm']
    if offtarget_calc_method not in possible_offtarget_calc_method:
        raise Exception('Unsupported offtarget calc method: ',offtarget_calc_method,'\nPlease select from the following: ',str(possible_offtarget_calc_method))
    possible_probe_selection_method = ['Greedy','Dynamic']
    if probe_selection_method not in possible_probe_selection_method:
        raise Exception('Unsupported probe selection method: ',possible_probe_selection_method,'\n Please select from the following: ', str(possible_probe_selection_method))
    
    # Print out config parameters 
    print('\tNumber of probes per target: ', str(num_probes))
    print('\tProbe length: ' + str(probe_len))
    print('\tDistance between probes: ' + str(probe_dist))
    print('\tMinimum probe Tm: ' + str(min_tm))
    print('\tMaximum probe hairpin Tm: ' + str(hairpin_max_tm))
    print('\tOfftarget calculation method: ' + str(offtarget_calc_method))
    print('\tProbe selection method: ' + str(probe_selection_method))
    print('\tFilter out repeat nucleotides?: ' + str(filter_repeats))
    
    # Set the transcriptome and database files using
    # the architecture of how the program runs
    transcriptome_file = 'species.'+species+'/transcriptome.fa'
    transcriptome_db = 'species.'+species+'/blast_db'

    # Make the temp directory if not already made
    if not os.path.exists('./temp/'):
        os.makedirs('./temp/')

    # Make the output folder if not already made
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    # Load transcriptome file into a dictionary with transcript_id as the key
    # Only loads if it wasn't passed into the main function by the GUI
    if transcriptome_dict == {}:
        transcriptome_dict, gene_info_array, rRNA_list = utils.load_transcriptome_dict(transcriptome_file)

    # Loads transcript targets if not loaded into the main function by the GUI
    # Does this using complex utils function
    if transcript_targets == {}:
        transcript_targets = utils.input_to_transcript_dictionary(transcriptome_dict, species,
                                   args.transcripts_file, args.gene_symbol_file, args.gene_symbols, args.ensembl_gene_ids, args.ensembl_transcript_ids)


    ##############################################################################
    #       Running main loop on each target to make probes                      #
    ##############################################################################

    # Sets variable to keep track of how many loops completed successfully
    successfully_made = 0
    
    # For loop which makes probes for all the inputted targets
    for target_transcript_id in transcript_targets:
        
        # Pulling gene ID and gene symbol from names in transcriptome file
        for gene_row in gene_info_array:
            if gene_row[1] == target_transcript_id:
                target_gene_id = gene_row[2]
                target_gene_symbol = gene_row[0]
        
        # Checking if ENSEMBL ID found, which will be used to exclude on target matches:
        try:
            target_gene_symbol
        except NameError:
            target_gene_id = 'unknown_id'
            target_gene_symbol = 'unknown_id'
            print("\nWarning: ENSEMBL ID was not found in input. Will proceed, but it may penalize the probes for binding its intended target.")

        print("Processing transcript " +
            target_gene_symbol + "_"+target_transcript_id + " ...\n")

        # Makes all probe possibilities and filters out ones with low Tm, high hairpin Tm, and repetitive sequences
        filtered_probe_tuples, probe_seq_dict = utils.create_and_filter_probe_possibilites(
            transcript_targets[target_transcript_id], probe_len, min_tm, hairpin_max_tm, num_probes, filter_repeats)

        # Write the filtered probe sequences to a fasta file for alignment via blast
        temp_probe_fasta_filename = './temp/' + target_gene_id + '_' + target_transcript_id + ".fa"
        utils.save_probe_file(filtered_probe_tuples, temp_probe_fasta_filename)

        # Alignment of the filtered probes by blast
        aln_data, blast_run = utils.align_blast(temp_probe_fasta_filename, transcriptome_db, probe_len)

        print("Alignments: "+str(len(aln_data)))
        print('Blast input args: \n'+' '.join(blast_run.args)+'\n')

        # Combining off target penalties from blast alignment
        # Using the blast alignment scores (extremely fast, but questionably accurate)
        if offtarget_calc_method == "Alignment_score":
            offtarget_scores_list = utils.compute_offtarget_scores(
                aln_data, target_transcript_id, rRNA_list)
        # Using Tm calculated by TmNN using primer3 calcHeterodimer (slow, but accurate and now MULTITHREADED)
        elif offtarget_calc_method == "Tm":
            offtarget_scores_list = utils.multithread_offtarget_tm(
                aln_data, target_transcript_id, rRNA_list)

        # Makes an output file of all the probes off target scores
        utils.save_scores_file(output_folder,target_transcript_id, target_gene_symbol, offtarget_scores_list, probe_seq_dict)
        # Select the final set of probes

        # Using the greedy method
        if probe_selection_method == "Greedy":
            print('Probe set selection being done by Greedy selection (picking the best probe each time)')
            final_probes = utils.select_nonoverlapping_probes(
                offtarget_scores_list, probe_seq_dict, probe_len, probe_dist, num_probes, min_tm)

        # Write the final list of probes to output file
        utils.save_output_file(final_probes, output_folder, target_gene_symbol,
                            target_gene_id, target_transcript_id, probe_seq_dict, output_format)
        
        # If it gets here then the probes were successfully made and the variable is incremented.
        successfully_made += 1
    # Deletes the temp folder 
    shutil.rmtree('./temp/')
    # Returns True or False for if the number of targets is equal to the number of probe files made.
    # This is used by the GUI to let the user know if it worked.
    return(len(transcript_targets) == successfully_made)


# Now that the main function is defined, if this script is the one being called as in `python COD_FISH.py`
# then it will run the main function
if __name__ == "__main__":
    main()
    
    # in any python script, as it's being run there are
    # variables that are automatically set, for example,
    # if this is the script that was run by `python COD_FISH.py`
    # in the command line then it will define __name__ as "__main__".
    # They follow the pattern of __something__