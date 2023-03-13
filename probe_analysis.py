from pathlib import Path
import sys
import os
import subprocess
import pprint
import matplotlib.pyplot as plt
import numpy as np

def save_probe_file(filtered_probe_array, temp_probe_fasta_filename):

    # Makes the probe file for alignment with blast
    probe_file = open(temp_probe_fasta_filename, "w")
    # Then for each probe outputs it in fasta format
    for probe_num, probe_seq in enumerate(filtered_probe_array):
        probe_name = 'probe_'+str(probe_num)
        probe_seq = probe_seq
        probe_file.write(">" + probe_name + "\n")
        probe_file.write(probe_seq + "\n")
    probe_file.close()
    return(temp_probe_fasta_filename)

def align_blast(blast_path, probe_fasta_filename, transcriptome_db, probe_len):
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
    
    blast_run = subprocess.run([blast_path,
                            '-task', task,
                            '-db', transcriptome_db,
                            '-query', probe_fasta_filename, 
                            '-num_threads', str(os.cpu_count()), 
                            '-strand', 'minus',
                            
                            '-outfmt','6 qacc sacc qseq sseq score'],
                            capture_output=True, 
                            creationflags = subprocess.CREATE_NO_WINDOW)

    raw_data = blast_run.stdout.decode("utf-8").split('\n')    # Converts to string then to a list where each element is a line
    raw_data.remove('')
    alignment_data = []
    for line in raw_data:
        alignment_data.append(line.split('\t'))               # Splits each alignment line at the tabs then adds it to the alignment_data list
    return(alignment_data, blast_run)

oligo_miner_folder = 'D:/Research/Thesis/test_COD-FISH/test_targets/oligominer_output/final/'

raj_folder = 'D:/Research/Thesis/test_COD-FISH/test_targets/raj_output/formatted/'

probedealer_folder = 'D:/Research/Thesis/test_COD-FISH/test_targets/probedealer_output/'

cod_folder = 'D:/Research/Thesis/test_COD-FISH/test_targets/cod_output/'

target_folder = 'D:/Research/Thesis/test_COD-FISH/test_targets/'


om_dict = {}

raj_dict = {}

pd_dict = {}

cod_dict = {}

source_dir = Path(oligo_miner_folder)
files = source_dir.glob('*.bed')
for file in files:
    
    target = os.path.basename(file).split('_probes')[0]

    om_dict[target] = []

    with open(file) as bed_file:
        for line in bed_file:
            om_dict[target].append(line.split('\t')[3])

for target in om_dict:
    if len(om_dict[target]) > 50:
        om_dict[target] = om_dict[target][0:50]

print(om_dict.keys)

source_dir = Path(raj_folder)
files = source_dir.glob('*.fa')
for file in files:
    
    target = os.path.basename(file).split('.')[0]

    raj_dict[target] = []

    with open(file) as fasta_file:
        for line in fasta_file:
            if '>' in line:
                continue
            else:
                raj_dict[target].append(line.strip())

print(raj_dict.keys())

source_dir = Path(probedealer_folder)
files = source_dir.glob('*.fasta')
for file in files:

    target = os.path.basename(file).split('.')[0]

    pd_dict[target] = []

    with open(file) as fasta_file:
        for line in fasta_file:
            if '>' in line:
                continue
            elif len(line) > 10:
                pd_dict[target].append(line.strip())

for target in pd_dict:
    if len(pd_dict[target]) > 50:
        pd_dict[target] = pd_dict[target][0:50]

print(pd_dict.keys())

source_dir = Path(cod_folder)
files = source_dir.glob('*.fa')
for file in files:

    target = os.path.basename(file).split('.')[0]

    target = target.split('_id_')[-1]

    cod_dict[target] = []

    with open(file) as fasta_file:
        for line in fasta_file:
            if '>' in line:
                continue
            elif len(line) > 10:
                cod_dict[target].append(line.strip())
            

print(cod_dict.keys())

#######################
# Number of probes per target

def len_dict(dictionary):
    len_dict = {}
    for target in dictionary:
        len_dict[target] = len(dictionary[target])
    return(len_dict)

def save_len_dict(len_dict, filename):
    with open(filename, 'w') as csv:
        for target in len_dict:
            csv.write(target+','+str(len_dict[target])+'\n')


om_len_dict = len_dict(om_dict)
pd_len_dict = len_dict(pd_dict)
raj_len_dict = len_dict(raj_dict)
cod_len_dict = len_dict(cod_dict)

save_len_dict(om_len_dict, 'om_len.csv')
save_len_dict(pd_len_dict, 'pd_len.csv')
save_len_dict(raj_len_dict, 'raj_len.csv')
save_len_dict(cod_len_dict, 'cod_len.csv')


#######################
# Plot lengths

# get averages:

def avg_lens(dictionary, seq_types):
    lens_by_type = {}
    for type in seq_types:
        lens_by_type[type] = []
        for target in dictionary:
            if type in target:
                lens_by_type[type].append(dictionary[target])
    avg_lens = {}
    for type in lens_by_type:
        try:
            avg_lens[type] = sum(lens_by_type[type])/len(lens_by_type[type])
        except ZeroDivisionError:
            print('zero len')
    return(avg_lens)

seq_types = ['Normal', 'ShortSeq', 'LongSeq', 'LowGC', 'HighGC']

om_avg_len = avg_lens(om_len_dict, seq_types)
pd_avg_len = avg_lens(pd_len_dict, seq_types)
raj_avg_len = avg_lens(raj_len_dict, seq_types)
cod_avg_len = avg_lens(cod_len_dict, seq_types)

om_avg_len_array = []
for type in seq_types:
    om_avg_len_array.append(om_avg_len[type])

pd_avg_len_array = []
for type in seq_types:
    pd_avg_len_array.append(pd_avg_len[type])

raj_avg_len_array = []
for type in seq_types:
    raj_avg_len_array.append(raj_avg_len[type])

cod_avg_len_array = []
for type in seq_types:
    cod_avg_len_array.append(cod_avg_len[type])

print(cod_avg_len)



x = np.arange(len(seq_types))
width = 0.2
#for type in seq_types:
plt.bar(x-3*width/2, om_avg_len_array,width, label='OligoMiner')
plt.bar(x-width/2, pd_avg_len_array,width, label='ProbeDealer')
plt.bar(x+width/2, raj_avg_len_array,width, label='RajLab')
plt.bar(x+3*width/2, cod_avg_len_array,width, label='COD-FISH')

plt.xticks(x, seq_types)
plt.ylabel("Average Number of Probes Generated")
plt.xlabel("RNA Target Category")

plt.ylim(0,70)
plt.legend(loc='upper left')
plt.show()

#######################
# Alignments per target

def align_sets(probe_dict):

    probe_set_scores = {}

    for target in probe_dict:
        # save probes
        temp_probe_fasta_filename = save_probe_file(probe_dict[target], "temp.fa")
        # align them
        alignment_data, blast_run = align_blast('blastn', temp_probe_fasta_filename, '../current/species.homo_sapiens/blast_db', 20)
        # score alignments
        probe_set_scores[target] = 0

        for aln in alignment_data:
            probe_set_scores[target] += int(aln[4])
    
    return(probe_set_scores)
'''
om_score_dict = align_sets(om_dict)
pd_score_dict = align_sets(pd_dict)
raj_score_dict = align_sets(raj_dict)
cod_score_dict = align_sets(cod_dict)

print(om_score_dict)
print(pd_score_dict)
print(raj_score_dict)
print(cod_score_dict)'''

om_score_dict = {'HighGC_0_Round1': 39523, 'HighGC_0_Round2': 13647, 'HighGC_0_Round3': 6953, 'HighGC_0_Round4': 9485, 'HighGC_0_Round5': 11094, 'HighGC_1_Round1': 14141, 'HighGC_1_Round2': 9203, 'HighGC_1_Round3': 18039, 'HighGC_1_Round4': 104, 'HighGC_1_Round5': 25947, 'HighGC_2_Round1': 8964, 'HighGC_2_Round2': 1479, 'HighGC_2_Round3': 14566, 'HighGC_2_Round4': 1002, 'HighGC_2_Round5': 23889, 'HighGC_3_Round1': 17628, 'HighGC_3_Round2': 25710, 'HighGC_3_Round3': 0, 'HighGC_3_Round4': 36773, 'HighGC_3_Round5': 28011, 'HighGC_4_Round1': 19950, 'HighGC_4_Round2': 23967, 'HighGC_4_Round3': 50701, 'HighGC_4_Round4': 31086, 'HighGC_4_Round5': 18856, 'LongSeq_0_Round1': 20177, 'LongSeq_0_Round2': 28603, 'LongSeq_0_Round3': 13451, 'LongSeq_0_Round4': 17963, 'LongSeq_0_Round5': 36543, 'LongSeq_1_Round1': 16627, 'LongSeq_1_Round2': 15832, 'LongSeq_1_Round3': 25111, 'LongSeq_1_Round4': 22091, 'LongSeq_1_Round5': 19697, 'LongSeq_2_Round1': 27138, 'LongSeq_2_Round2': 31762, 'LongSeq_2_Round3': 31945, 'LongSeq_2_Round4': 15224, 'LongSeq_2_Round5': 12933, 'LongSeq_3_Round1': 24855, 'LongSeq_3_Round2': 10022, 'LongSeq_3_Round3': 42808, 'LongSeq_3_Round4': 38164, 'LongSeq_3_Round5': 24154, 'LongSeq_4_Round1': 13131, 'LongSeq_4_Round2': 12376, 'LongSeq_4_Round3': 27601, 'LongSeq_4_Round4': 14432, 'LongSeq_4_Round5': 16891, 'LowGC_0_Round1': 26473, 'LowGC_0_Round2': 20570, 'LowGC_0_Round3': 15904, 'LowGC_0_Round4': 23096, 'LowGC_0_Round5': 39809, 'LowGC_1_Round1': 19122, 'LowGC_1_Round2': 15077, 'LowGC_1_Round3': 17335, 'LowGC_1_Round4': 17189, 'LowGC_1_Round5': 15222, 'LowGC_2_Round1': 11166, 'LowGC_2_Round2': 14572, 'LowGC_2_Round3': 13587, 'LowGC_2_Round4': 16616, 'LowGC_2_Round5': 17496, 'LowGC_3_Round1': 29327, 'LowGC_3_Round2': 19776, 'LowGC_3_Round3': 17129, 'LowGC_3_Round4': 17058, 'LowGC_3_Round5': 17638, 'LowGC_4_Round1': 14195, 'LowGC_4_Round2': 15375, 'LowGC_4_Round3': 22932, 'LowGC_4_Round4': 14739, 'LowGC_4_Round5': 18537, 'Normal_0': 0, 'Normal_1': 28252, 'Normal_10': 21258, 'Normal_11': 20755, 'Normal_12': 14966, 'Normal_14': 21549, 'Normal_15': 14567, 'Normal_16': 22540, 'Normal_18': 11941, 'Normal_19': 16240, 'Normal_2': 12642, 'Normal_20': 14946, 'Normal_21': 8328, 'Normal_22': 40366, 'Normal_23': 27002, 'Normal_24': 24311, 'Normal_25': 15072, 'Normal_26': 14562, 'Normal_27': 49209, 'Normal_3': 11207, 'Normal_4': 18350, 'Normal_5': 31549, 'Normal_7': 24088, 'Normal_8': 19865, 'Normal_9': 16023, 'ShortSeq_0_Round1': 25951, 'ShortSeq_0_Round2': 29, 'ShortSeq_0_Round3': 25923, 'ShortSeq_0_Round4': 14967, 'ShortSeq_0_Round5': 8272, 'ShortSeq_1_Round1': 0, 'ShortSeq_1_Round2': 17959, 'ShortSeq_1_Round3': 10673, 'ShortSeq_1_Round4': 17029, 'ShortSeq_1_Round5': 10755, 'ShortSeq_2_Round1': 16673, 'ShortSeq_2_Round2': 15945, 'ShortSeq_2_Round3': 5828, 'ShortSeq_2_Round4': 26770, 'ShortSeq_2_Round5': 306, 'ShortSeq_3_Round1': 1952, 'ShortSeq_3_Round2': 31175, 'ShortSeq_3_Round3': 14727, 'ShortSeq_3_Round4': 45662, 'ShortSeq_3_Round5': 21937, 'ShortSeq_4_Round1': 17459, 'ShortSeq_4_Round2': 14569, 'ShortSeq_4_Round3': 17611, 'ShortSeq_4_Round4': 16783, 'ShortSeq_4_Round5': 16197}
pd_score_dict = {'ShortSeq_1_Round5': 13058, 'HighGC_0_Round2': 13535, 'HighGC_0_Round4': 5322, 'HighGC_1_Round2': 38477, 'HighGC_1_Round4': 6808, 'LongSeq_0_Round2': 16177, 'LongSeq_1_Round4': 17964, 'LowGC_0_Round4': 16799, 'LowGC_1_Round2': 17973, 'LowGC_1_Round4': 18403, 'Normal_1': 24218, 'Normal_3': 10947, 'Normal_7': 17695, 'Normal_19': 17097, 'Normal_22': 20348, 'Normal_24': 15625, 'Normal_26': 14255, 'ShortSeq_2_Round1': 11184, 'ShortSeq_2_Round5': 12615, 'ShortSeq_3_Round1': 7894, 'ShortSeq_3_Round3': 6143, 'ShortSeq_4_Round1': 21522, 'ShortSeq_4_Round3': 25052, 'HighGC_0_Round1': 9058, 'HighGC_0_Round3': 18515, 'HighGC_0_Round5': 9923, 'HighGC_1_Round1': 10794, 'HighGC_1_Round3': 12471, 'HighGC_1_Round5': 4536, 'HighGC_2_Round5': 14700, 'HighGC_3_Round1': 14414, 'HighGC_4_Round1': 11419, 'LongSeq_0_Round1': 25048, 'LongSeq_0_Round3': 18502, 'LongSeq_0_Round4': 13334, 'LongSeq_0_Round5': 14571, 'LongSeq_1_Round1': 12095, 'LongSeq_1_Round2': 14341, 'LongSeq_1_Round3': 13780, 'LongSeq_1_Round5': 59394, 'LongSeq_2_Round5': 13774, 'LongSeq_3_Round1': 12600, 'LongSeq_3_Round3': 16449, 'LongSeq_4_Round1': 12842, 'LongSeq_4_Round3': 81099, 'LowGC_0_Round1': 15872, 'LowGC_0_Round2': 25501, 'LowGC_0_Round3': 11946, 'LowGC_0_Round5': 12978, 'LowGC_1_Round1': 51241, 'LowGC_1_Round3': 16737, 'LowGC_1_Round5': 16451, 'LowGC_2_Round1': 17586, 'LowGC_4_Round5': 11215, 'Normal_0': 9058, 'Normal_2': 16230, 'Normal_4': 34736, 'Normal_5': 17763, 'Normal_8': 50656, 'Normal_12': 17240, 'Normal_18': 14159, 'Normal_20': 11676, 'Normal_21': 10040, 'Normal_23': 18222, 'Normal_25': 29504, 'Normal_27': 10650, 'ShortSeq_0_Round2': 18035, 'ShortSeq_1_Round4': 66221, 'ShortSeq_2_Round2': 34819, 'ShortSeq_2_Round3': 77258, 'ShortSeq_2_Round4': 27029, 'ShortSeq_3_Round2': 7421, 'ShortSeq_3_Round5': 10475, 'ShortSeq_4_Round2': 8809, 'ShortSeq_4_Round4': 16085, 'ShortSeq_4_Round5': 11643, 'HighGC_2_Round1': 4608, 'HighGC_2_Round2': 7063, 'HighGC_2_Round3': 7253, 'HighGC_2_Round4': 12992, 'HighGC_3_Round2': 6471, 'HighGC_3_Round3': 5395, 'HighGC_3_Round4': 21552, 'HighGC_3_Round5': 8081, 'HighGC_4_Round2': 8015, 'HighGC_4_Round3': 23040, 'HighGC_4_Round4': 7361, 'HighGC_4_Round5': 8081, 'LongSeq_2_Round1': 15736, 'LongSeq_2_Round2': 15261, 'LongSeq_2_Round3': 15548, 'LongSeq_2_Round4': 19908, 'LongSeq_3_Round2': 16619, 'LongSeq_3_Round4': 19850, 'LongSeq_3_Round5': 14570, 'LongSeq_4_Round2': 12821, 'LongSeq_4_Round5': 17086, 'LowGC_2_Round2': 16369, 'LowGC_2_Round3': 14419, 'LowGC_2_Round4': 14515, 'LowGC_2_Round5': 11997, 'LowGC_3_Round1': 14523, 'LowGC_3_Round2': 13857, 'LowGC_3_Round3': 13420, 'LowGC_3_Round4': 15406, 'LowGC_3_Round5': 16803, 'LowGC_4_Round1': 26103, 'LowGC_4_Round2': 14174, 'LowGC_4_Round3': 15434, 'LowGC_4_Round4': 11476, 'Normal_9': 18973, 'Normal_10': 13726, 'Normal_11': 13411, 'Normal_14': 14967, 'Normal_15': 19858, 'Normal_16': 38917, 'ShortSeq_0_Round1': 11582, 'ShortSeq_0_Round3': 10440, 'ShortSeq_0_Round4': 21187, 'ShortSeq_0_Round5': 12989, 'ShortSeq_1_Round1': 5148, 'ShortSeq_1_Round2': 17986, 'ShortSeq_1_Round3': 13233}
raj_score_dict = {'HighGC_0_Round1': 29951, 'HighGC_0_Round2': 33211, 'HighGC_0_Round3': 44214, 'HighGC_0_Round4': 40786, 'HighGC_0_Round5': 23801, 'HighGC_1_Round1': 39802, 'HighGC_1_Round2': 61858, 'HighGC_1_Round3': 35572, 'HighGC_1_Round4': 210010, 'HighGC_1_Round5': 34285, 'HighGC_2_Round1': 36974, 'HighGC_2_Round2': 63085, 'HighGC_2_Round3': 36232, 'HighGC_2_Round4': 64338, 'HighGC_2_Round5': 31849, 'HighGC_3_Round1': 27016, 'HighGC_3_Round2': 27462, 'HighGC_3_Round3': 33112, 'HighGC_3_Round4': 50612, 'HighGC_3_Round5': 27216, 'HighGC_4_Round1': 18308, 'HighGC_4_Round2': 35601, 'HighGC_4_Round3': 59123, 'HighGC_4_Round4': 30472, 'HighGC_4_Round5': 24430, 'LongSeq_0_Round1': 51978, 'LongSeq_0_Round2': 26373, 'LongSeq_0_Round3': 39229, 'LongSeq_0_Round4': 21037, 'LongSeq_0_Round5': 32758, 'LongSeq_1_Round1': 39477, 'LongSeq_1_Round2': 56250, 'LongSeq_1_Round3': 23901, 'LongSeq_1_Round4': 19978, 'LongSeq_1_Round5': 59138, 'LongSeq_2_Round1': 42205, 'LongSeq_2_Round2': 32311, 'LongSeq_2_Round3': 32704, 'LongSeq_2_Round4': 21141, 'LongSeq_2_Round5': 13673, 'LongSeq_3_Round1': 25030, 'LongSeq_3_Round2': 65217, 'LongSeq_3_Round3': 39657, 'LongSeq_3_Round4': 30216, 'LongSeq_3_Round5': 27660, 'LongSeq_4_Round1': 29520, 'LongSeq_4_Round2': 16017, 'LongSeq_4_Round3': 63956, 'LongSeq_4_Round4': 39267, 'LongSeq_4_Round5': 49519, 'LowGC_0_Round1': 32433, 'LowGC_0_Round2': 25536, 'LowGC_0_Round3': 19907, 'LowGC_0_Round4': 27078, 'LowGC_0_Round5': 35392, 'LowGC_1_Round1': 42518, 'LowGC_1_Round2': 283915, 'LowGC_1_Round3': 18947, 'LowGC_1_Round4': 73076, 'LowGC_1_Round5': 14609, 'LowGC_2_Round1': 24338, 'LowGC_2_Round2': 23150, 'LowGC_2_Round3': 18221, 'LowGC_2_Round4': 28994, 'LowGC_2_Round5': 27000, 'LowGC_3_Round1': 29882, 'LowGC_3_Round2': 13741, 'LowGC_3_Round3': 32785, 'LowGC_3_Round4': 19168, 'LowGC_3_Round5': 18856, 'LowGC_4_Round1': 27146, 'LowGC_4_Round2': 16606, 'LowGC_4_Round3': 17756, 'LowGC_4_Round4': 15387, 'LowGC_4_Round5': 19822, 'Normal_0': 118897, 'Normal_1': 46953, 'Normal_10': 20491, 'Normal_11': 18570, 'Normal_12': 21463, 'Normal_14': 22095, 'Normal_15': 34341, 'Normal_16': 69014, 'Normal_18': 24771, 'Normal_19': 21333, 'Normal_2': 19981, 'Normal_20': 22748, 'Normal_21': 32436, 'Normal_22': 38414, 'Normal_23': 24750, 'Normal_24': 19018, 'Normal_25': 53780, 'Normal_26': 34733, 'Normal_27': 47380, 'Normal_3': 42822, 'Normal_4': 47259, 'Normal_5': 33517, 'Normal_7': 19660, 'Normal_8': 60177, 'Normal_9': 37400, 'ShortSeq_0_Round1': 33991, 'ShortSeq_0_Round2': 68184, 'ShortSeq_0_Round3': 23708, 'ShortSeq_0_Round4': 24723, 'ShortSeq_0_Round5': 23433, 'ShortSeq_1_Round1': 25593, 'ShortSeq_1_Round2': 24395, 'ShortSeq_1_Round3': 13223, 'ShortSeq_1_Round4': 94863, 'ShortSeq_1_Round5': 25621, 'ShortSeq_2_Round1': 41610, 'ShortSeq_2_Round2': 87995, 'ShortSeq_2_Round3': 77912, 'ShortSeq_2_Round4': 51637, 'ShortSeq_2_Round5': 27926, 'ShortSeq_3_Round1': 120696, 'ShortSeq_3_Round2': 31617, 'ShortSeq_3_Round3': 14901, 'ShortSeq_3_Round4': 37016, 'ShortSeq_3_Round5': 17436, 'ShortSeq_4_Round1': 33946, 'ShortSeq_4_Round2': 23488, 'ShortSeq_4_Round3': 141048, 'ShortSeq_4_Round4': 35588, 'ShortSeq_4_Round5': 30813}
cod_score_dict = {'HighGC_0_Round1': 10628, 'HighGC_0_Round2': 11841, 'HighGC_0_Round3': 10038, 'HighGC_0_Round4': 17763, 'HighGC_0_Round5': 11801, 'HighGC_1_Round1': 2161, 'HighGC_1_Round2': 24313, 'HighGC_1_Round3': 12047, 'HighGC_1_Round4': 140189, 'HighGC_1_Round5': 10331, 'HighGC_2_Round1': 13983, 'HighGC_2_Round2': 37643, 'HighGC_2_Round3': 10636, 'HighGC_2_Round4': 41081, 'HighGC_2_Round5': 6760, 'HighGC_3_Round1': 7270, 'HighGC_3_Round2': 4719, 'HighGC_3_Round3': 13444, 'HighGC_3_Round4': 12317, 'HighGC_3_Round5': 8200, 'HighGC_4_Round1': 5431, 'HighGC_4_Round2': 10158, 'HighGC_4_Round3': 30528, 'HighGC_4_Round4': 8923, 'HighGC_4_Round5': 2481, 'LongSeq_0_Round1': 4907, 'LongSeq_0_Round2': 12800, 'LongSeq_0_Round3': 6571, 'LongSeq_0_Round4': 5924, 'LongSeq_0_Round5': 12914, 'LongSeq_1_Round1': 15830, 'LongSeq_1_Round2': 4824, 'LongSeq_1_Round3': 1903, 'LongSeq_1_Round4': 2942, 'LongSeq_1_Round5': 3792, 'LongSeq_2_Round1': 24200, 'LongSeq_2_Round2': 5942, 'LongSeq_2_Round3': 12213, 'LongSeq_2_Round4': 2073, 'LongSeq_2_Round5': 2361, 'LongSeq_3_Round1': 11935, 'LongSeq_3_Round2': 5407, 'LongSeq_3_Round3': 9466, 'LongSeq_3_Round4': 6494, 'LongSeq_3_Round5': 4780, 'LongSeq_4_Round1': 5306, 'LongSeq_4_Round2': 2762, 'LongSeq_4_Round3': 1931, 'LongSeq_4_Round4': 1924, 'LongSeq_4_Round5': 36804, 'LowGC_0_Round1': 4247, 'LowGC_0_Round2': 1530, 'LowGC_0_Round3': 2263, 'LowGC_0_Round4': 2724, 'LowGC_0_Round5': 13712, 'LowGC_1_Round1': 7610, 'LowGC_1_Round2': 3806, 'LowGC_1_Round3': 2226, 'LowGC_1_Round4': 42785, 'LowGC_1_Round5': 5044, 'LowGC_2_Round1': 8656, 'LowGC_2_Round2': 3919, 'LowGC_2_Round3': 2305, 'LowGC_2_Round4': 8582, 'LowGC_2_Round5': 7275, 'LowGC_3_Round1': 9438, 'LowGC_3_Round2': 1430, 'LowGC_3_Round3': 8008, 'LowGC_3_Round4': 10133, 'LowGC_3_Round5': 6460, 'LowGC_4_Round1': 9326, 'LowGC_4_Round2': 4891, 'LowGC_4_Round3': 4966, 'LowGC_4_Round4': 3151, 'LowGC_4_Round5': 7072, 'Normal_0': 52682, 'Normal_1': 13599, 'Normal_10': 9512, 'Normal_11': 4305, 'Normal_12': 7767, 'Normal_14': 3574, 'Normal_15': 16794, 'Normal_16': 5470, 'Normal_18': 8647, 'Normal_19': 6600, 'Normal_2': 3922, 'Normal_20': 4145, 'Normal_21': 14334, 'Normal_22': 11169, 'Normal_23': 6991, 'Normal_24': 5127, 'Normal_25': 12298, 'Normal_26': 9658, 'Normal_27': 23785, 'Normal_3': 24466, 'Normal_4': 2671, 'Normal_5': 9654, 'Normal_7': 6100, 'Normal_8': 4250, 'Normal_9': 23645, 'ShortSeq_0_Round1': 19459, 'ShortSeq_0_Round2': 53881, 'ShortSeq_0_Round3': 11155, 'ShortSeq_0_Round4': 4251, 'ShortSeq_0_Round5': 10247, 'ShortSeq_1_Round1': 15087, 'ShortSeq_1_Round2': 12323, 'ShortSeq_1_Round3': 3566, 'ShortSeq_1_Round4': 7326, 'ShortSeq_1_Round5': 11549, 'ShortSeq_2_Round1': 27151, 'ShortSeq_2_Round2': 36359, 'ShortSeq_2_Round3': 26128, 'ShortSeq_2_Round4': 15507, 'ShortSeq_2_Round5': 18044, 'ShortSeq_3_Round1': 65622, 'ShortSeq_3_Round2': 12629, 'ShortSeq_3_Round3': 5887, 'ShortSeq_3_Round4': 26091, 'ShortSeq_3_Round5': 9936, 'ShortSeq_4_Round1': 6361, 'ShortSeq_4_Round2': 15475, 'ShortSeq_4_Round3': 100305, 'ShortSeq_4_Round4': 13046, 'ShortSeq_4_Round5': 7870}


def normalize(score_dict,len_dict):
    norm_score_dict = {}
    for target in score_dict:
        try:
            norm_score_dict[target] = score_dict[target] / len_dict[target]
        except ZeroDivisionError:
            print('zero div')

    return(norm_score_dict)


avg_align_om = avg_lens(om_score_dict, seq_types)
avg_align_pd = avg_lens(pd_score_dict, seq_types)
avg_align_raj = avg_lens(raj_score_dict, seq_types)
avg_align_cod = avg_lens(cod_score_dict, seq_types)




om_avg_score_array = []
for type in seq_types:
    om_avg_score_array.append(avg_align_om[type])

pd_avg_score_array = []
for type in seq_types:
    pd_avg_score_array.append(avg_align_pd[type])

raj_avg_score_array = []
for type in seq_types:
    raj_avg_score_array.append(avg_align_raj[type])

cod_avg_score_array = []
for type in seq_types:
    cod_avg_score_array.append(avg_align_cod[type])


om_per_probe_score = [score / num for score,num in zip(om_avg_score_array, om_avg_len_array)]
pd_per_probe_score = [score / num for score,num in zip(pd_avg_score_array, pd_avg_len_array)]
raj_per_probe_score = [score / num for score,num in zip(raj_avg_score_array, raj_avg_len_array)]
cod_per_probe_score = [score / num for score,num in zip(cod_avg_score_array, cod_avg_len_array)]


x = np.arange(len(seq_types))
width = 0.2
#for type in seq_types:
'''plt.bar(x-3*width/2, om_avg_score_array,width, label='OligoMiner')
plt.bar(x-width/2, pd_avg_score_array,width, label='ProbeDealer')
plt.bar(x+width/2, raj_avg_score_array,width, label='RajLab')
plt.bar(x+3*width/2, cod_avg_score_array,width, label='COD-FISH')'''

plt.bar(x-3*width/2, om_per_probe_score,width, label='OligoMiner')
plt.bar(x-width/2, pd_per_probe_score,width, label='ProbeDealer')
plt.bar(x+width/2, raj_per_probe_score,width, label='RajLab')
plt.bar(x+3*width/2, cod_per_probe_score,width, label='COD-FISH')

plt.xticks(x, seq_types)
plt.ylabel("Average Cumulative Alignment Score /\n Average Number of Probes")
plt.xlabel("RNA Target Category")
plt.ylim(0,1100)
plt.legend(loc='upper left')
plt.show()