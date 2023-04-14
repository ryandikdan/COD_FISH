# This script will be used to make the functions for the graphical user interface (GUI) using tkinter
# This script will then be combined and packaged as an exe for windows or as a program for linux/macOS

# Used Pyinstaller to compile script to .exe in Powershell with the following command:
# `pyinstaller.exe -i .\codfish_icon.png --add-data='codfish_icon.png;.' .\gui_COD_FISH.py`
# -i for icon -w for windowed (no terminal) and -F for making it a single file executable
# This includes the icon file and uses it for the program that's made.
# ! This doesn't work completely and doesn't add the primer3 or ensembl_rest folders
# completely to the dist folder, but upon adding them manually it worked fine.
# Comtemplating allowing the terminal to be created so the process can be monitored.

# For mac or linux have to use this command:
# `pyinstaller -i codfish_icon.png --add-data='codfish_icon.png:.' gui_COD_FISH.py`

# Trying to import in only the necessary parts of Tkinter so that pyinstaller executable is not too laggy

from tkinter import ttk, Tk, Entry, PhotoImage, filedialog, messagebox, END, INSERT

import ensembl_rest
import multiprocessing as mp
import utils_COD_FISH as utils
import COD_FISH
import subprocess
import os
import platform


# This awesome script was found online which makes the finding of species via autocomplete possible.
# All credit to these wonderful people for their fantastic code!:  
# Mitja Martini and Russell Adams
# Created by Mitja Martini on 2008-11-29.
# Updated by Russell Adams, 2011/01/24 to support Python 3 and Combobox.
#    Licensed same as original (not specified?), or public domain, whichever is less restrictive.

class AutocompleteEntry(Entry):
    """
    Subclass of Tkinter.Entry that features autocompletion.

    To enable autocompletion use set_completion_list(list) to define
    a list of possible strings to hit.
    To cycle through hits use down and up arrow keys.
    """

    def set_completion_list(self, completion_list):
        self._completion_list = sorted(
            completion_list, key=str.lower)  # Work with a sorted list
        self._hits = []
        self._hit_index = 0
        self.position = 0
        self.bind('<KeyRelease>', self.handle_keyrelease)

    def autocomplete(self, delta=0):
        """autocomplete the Entry, delta may be 0/1/-1 to cycle through possible hits"""
        if delta:  # need to delete selection otherwise we would fix the current position
            self.delete(self.position, END)
        else:  # set position to end so selection starts where textentry ended
            self.position = len(self.get())
        # collect hits
        _hits = []
        for element in self._completion_list:
            if element.lower().startswith(self.get().lower()):  # Match case-insensitively
                _hits.append(element)
        # if we have a new hit list, keep this in mind
        if _hits != self._hits:
            self._hit_index = 0
            self._hits = _hits
        # only allow cycling if we are in a known hit list
        if _hits == self._hits and self._hits:
            self._hit_index = (self._hit_index + delta) % len(self._hits)
        # now finally perform the auto completion
        if self._hits:
            self.delete(0, END)
            self.insert(0, self._hits[self._hit_index])
            self.select_range(self.position, END)

    def handle_keyrelease(self, event):
        # I may have broken the unicode...
        tkinter_umlauts = ['odiaeresis', 'adiaeresis', 'udiaeresis',
                           'Odiaeresis', 'Adiaeresis', 'Udiaeresis', 'ssharp']
        """event handler for the keyrelease event on this widget"""
        if event.keysym == "BackSpace":
            self.delete(self.index(INSERT), END)
            self.position = self.index(END)
        if event.keysym == "Left":
            if self.position < self.index(END):  # delete the selection
                self.delete(self.position, END)
            else:
                self.position = self.position-1  # delete one character
                self.delete(self.position, END)
        if event.keysym == "Right":
            self.position = self.index(END)  # go to end (no selection)
        if event.keysym == "Down":
            self.autocomplete(1)  # cycle to next hit
        if event.keysym == "Up":
            self.autocomplete(-1)  # cycle to previous hit
        if len(event.keysym) == 1 or event.keysym in tkinter_umlauts:
            self.autocomplete()

class AutocompleteCombobox(ttk.Combobox):

    def set_completion_list(self, completion_list):
        """Use our completion list as our drop down selection menu, arrows move through menu."""
        self._completion_list = sorted(
            completion_list, key=str.lower)  # Work with a sorted list
        self._hits = []
        self._hit_index = 0
        self.position = 0
        self.bind('<KeyRelease>', self.handle_keyrelease)
        self['values'] = self._completion_list  # Setup our popup menu

    def autocomplete(self, delta=0):
        """autocomplete the Combobox, delta may be 0/1/-1 to cycle through possible hits"""
        if delta:  # need to delete selection otherwise we would fix the current position
            self.delete(self.position, END)
        else:  # set position to end so selection starts where textentry ended
            self.position = len(self.get())
        # collect hits
        _hits = []
        for element in self._completion_list:
            if element.lower().startswith(self.get().lower()):  # Match case insensitively
                _hits.append(element)
        # if we have a new hit list, keep this in mind
        if _hits != self._hits:
            self._hit_index = 0
            self._hits = _hits
        # only allow cycling if we are in a known hit list
        if _hits == self._hits and self._hits:
            self._hit_index = (self._hit_index + delta) % len(self._hits)
        # now finally perform the auto completion
        if self._hits:
            self.delete(0, END)
            self.insert(0, self._hits[self._hit_index])
            self.select_range(self.position, END)

    def handle_keyrelease(self, event):
        """event handler for the keyrelease event on this widget"""
        if event.keysym == "BackSpace":
            self.delete(self.index(INSERT), END)
            self.position = self.index(END)
        if event.keysym == "Left":
            if self.position < self.index(END):  # delete the selection
                self.delete(self.position, END)
            else:
                self.position = self.position-1  # delete one character
                self.delete(self.position, END)
        if event.keysym == "Right":
            self.position = self.index(END)  # go to end (no selection)
        if len(event.keysym) == 1:
            self.autocomplete()
        # No need for up/down, we'll jump to the popup
        # list at the position of the autocompletion
        
# I'm not a huge fan of making the whole script a main function
# but in order to avoid any weird interactions between the 
# multithreading of multiprocessing opening a million tkinter windows 
# this should stop that.
def main():
    
    # Pull the species names in the format 'homo_sapiens' from ensembl REST API
    ensembl_species = [species_dict['name'] for species_dict in ensembl_rest.species()['species']]

    # Converts 'homo_sapiens' to 'Homo sapiens' which is just easier on the eyes
    # May use display names going forward 
    formatted_species = [ensembl_species[species_idx].replace(
        '_', ' ').capitalize() for species_idx in range(len(ensembl_species))]

    # Make the default config file if it's not present
    if not os.path.exists("config.py"):
        utils.make_config_file()


    ###################################################
    #     Homemade tkinter helper functions
    ###################################################

    def build_input(type, default_val, text, inside_frame, col, row, next_to_or_under, width, sticky='ew', values=[], padding=5, format = ''):
        # Function to make the entries easier to code.
        # All this does is take what each entry, combobox, button, etc needs and then
        # runs the appropriate tkinter commands, returning the label and entry variables.
        # You still have to remember and know the rows and columns for input, but that's not hard.

        # This code makes the label and entry either next to or under each other
        next_to = 0
        under = 0
        if next_to_or_under == 'next_to':
            next_to = 1
        if next_to_or_under == 'under':
            under = 1

        # For text entry
        if type == 'entry':
            label_var = ttk.Label(inside_frame, text=text)
            label_var.grid(column=col, row=row, padx=padding,
                        pady=padding, sticky=sticky)

            entry_var = ttk.Entry(inside_frame, width=width)
            entry_var.insert(0, default_val)
            entry_var.grid(column=col+next_to, row=row+under,
                        padx=padding, pady=padding, sticky=sticky)
            return(label_var, entry_var)

        # For combobox entry (drop down menu)
        # Default value picks that number element as default
        if type == 'combobox':
            label_var = ttk.Label(inside_frame, text=text)
            label_var.grid(column=col, row=row, padx=padding,
                        pady=padding, sticky=sticky)

            combobox_var = ttk.Combobox(inside_frame, width=width, values=values)
            combobox_var.current(default_val)
            combobox_var.grid(column=col+next_to, row=row+under,
                            padx=padding, pady=padding, sticky=sticky)
            return(label_var, combobox_var)

        # For autocombobox entry (drop down menu with autocomplete on type)
        # Default value picks that number element as default
        if type == 'autocombobox':
            label_var = ttk.Label(inside_frame, text=text)
            label_var.grid(column=col, row=row, sticky=sticky)

            autocombobox_var = AutocompleteCombobox(inside_frame, width=width)
            autocombobox_var.set_completion_list(values)
            autocombobox_var.focus_set()
            autocombobox_var.grid(column=col+next_to, row=row +
                                under, padx=padding, pady=padding, sticky=sticky)
            return(label_var, autocombobox_var)
        
        # For file select entry
        # Puts the button below the text, can move the +1 under the buttons.grid(row) to change that.
        if type == 'filedialog':
            label_var = ttk.Label(inside_frame, text=text)
            label_var.grid(column=col, row=row)
            
            entry_var = ttk.Entry(inside_frame, width=width)
            entry_var.insert(0, default_val)
            entry_var.grid(column=col+next_to, row=row+under,
                        padx=padding, pady=padding, sticky=sticky)
            
            filedialog_button = ttk.Button(inside_frame, text='Select file...', command=lambda: select_file(entry_var, format))
            filedialog_button.grid(column=col+next_to, row=row+under+1, sticky=sticky, padx=padding, pady=padding)
            
            return(label_var, entry_var, filedialog_button)

        else:
            raise(Exception('Invalid type'))


    def build_label_frame(parent_frame, col, row, text):
        frame = ttk.LabelFrame(parent_frame, text=text)
        frame.grid(column=col, row=row, padx=10, pady=10, sticky='ew')
        return(frame)


    def get_display_size():
        root = Tk()
        root.update_idletasks()
        root.attributes('-fullscreen', True)
        root.state('iconic')
        height = root.winfo_screenheight()
        width = root.winfo_screenwidth()
        root.destroy()
        return(height, width)


    ###################################################
    #       Button functions
    ###################################################

    def run_program(event=None):
        # Events that happen when the main button is clicked.
        
        # Species
        formatted_species = species_autocombobox.get()
        species = formatted_species.lower().replace(' ','_')
        # all_species already defined at top
        species_brief = ['human', 'mouse', 'zebrafish', 'rat',
                        'pig', 'chimpanzee', 'chicken', 'cow', 'fruitfly', 'celegans']
        brief_to_ensembl_dict = {'human': 'homo_sapiens', 'mouse': 'mus_musculus', 'zebrafish': 'danio_rerio', 'rat': 'rattus_norvegicus',
                                'pig': 'sus_scrofa', 'chimpanzee': 'pan_troglodytes', 'chicken': 'gallus_gallus', 'cow': 'bos_taurus', 'fruitfly': 'drosophila_melanogaster', 'celegans': 'caenorhabditis_elegans'}
        if species in species_brief:
            species = brief_to_ensembl_dict[species]
        if (species not in ensembl_species) and (species not in species_brief):
            raise Exception("Error\nThe species that you have entered,", species, "cannot be found. \
                            Try to use autocomplete and the dropdown menu to find your species of interest and make sure it's spelled correctly.")
        
        # Genes
        if transcript_file_var.get() != '':
            transcripts_file = transcript_file_var.get()
        else:
            transcripts_file = None
        if gene_symbol_file_var.get() != '':
            gene_symbol_file = gene_symbol_file_var.get()
        else:
            gene_symbol_file = None
            
        # Using these if statements to split up the output based on how the user input it
        # for gene symbols
        gene_symbols = []
        if gene_symbol_var.get() != '':
            if ', ' in gene_symbol_var.get():
                gene_symbols = gene_symbol_var.get().split(', ')
            elif ' ' in gene_symbol_var.get():
                gene_symbols = gene_symbol_var.get().split(' ')
            elif ',' in gene_symbol_var.get():
                gene_symbols = gene_symbol_var.get().split(',')
            else:
                gene_symbols.append(gene_symbol_var.get())
        else:
            gene_symbols = []
        # for ensembl gene IDs
        ensembl_gene_ids = []
        if gene_ID_var.get() != '':
            if ', ' in gene_ID_var.get():
                ensembl_gene_ids = gene_ID_var.get().split(', ')
            elif ' ' in gene_ID_var.get():
                ensembl_gene_ids = gene_ID_var.get().split(' ')
            elif ',' in gene_ID_var.get():
                ensembl_gene_ids = gene_ID_var.get().split(',')
            else:
                ensembl_gene_ids.append(gene_ID_var.get())
        else:
            ensembl_gene_ids = []
        # for ensembl transcript IDs
        ensembl_transcript_ids = []
        if transcript_ID_var.get() != '':
            if ', ' in transcript_ID_var.get():
                ensembl_transcript_ids = transcript_ID_var.get().split(', ')
            elif ' ' in transcript_ID_var.get():
                ensembl_transcript_ids = transcript_ID_var.get().split(' ')
            elif ',' in transcript_ID_var.get():
                ensembl_transcript_ids = transcript_ID_var.get().split(',')
            else:
                ensembl_transcript_ids.append(transcript_ID_var.get())
        else:
            ensembl_transcript_ids = []

        # Checks and makes species folder if needed
        utils.make_and_check_species_folder(species)
        # Loads the transcriptome.fa file for the species indicated
        transcriptome_dict, gene_info_dataframe, rRNA_list = utils.load_transcriptome_dict('./species.'+species+'/transcriptome.fa')
        # Loads the targets for input into COD_FISH.main
        transcript_targets = utils.input_to_transcript_dictionary(transcriptome_dict, species,
                                                                    transcripts_file, gene_symbol_file, gene_symbols, ensembl_gene_ids, ensembl_transcript_ids)
        # Runs the actual program, which returns result, a boolean that says if the probes were generated as intended.
        result = COD_FISH.main(transcript_targets, species, transcriptome_dict, gene_info_dataframe, rRNA_list)
        
        if result:
            messagebox.showinfo('Probes were generated','Click OK to open the output folder for the generated probe files')
            # If windows, opens windows explorer in the output folder
            if os.name =='nt':
                subprocess.run(['explorer', os.getcwd()+'\\output\\'])
            elif platform.system() == 'Darwin':
                subprocess.call(["open", "-R", os.getcwd()+'/output/'])
        else:
            messagebox.showerror('Something went wrong','Please check the input parameters and if the issue can\'t be resolved then please raise it in the Github.')
        # !!! Kills window to avoid error of adding already made targets to new queries
        # !!! To fix later
        root.destroy()
        
    def alter_config(event=None):
        # This is the action that happens when the advanced pane's button is clicked
        
        # Removes config.py if it already exists
        if os.path.exists("config.py"):
            os.remove("config.py")
            
        # Sets the variables equal to what's in the GUI
        num_probes = probe_num_entry.get()
        probe_len = probe_len_entry.get()
        probe_dist = probe_dist_entry.get()
        min_tm = probe_min_tm_entry.get()
        hairpin_max_tm = probe_max_hairpin_tm_entry.get()
        offtarget_calc_method = off_target_method_entry.get()
        probe_selection_method = probe_selection_entry.get()
        filter_repeats = filter_repeats_entry.get()
        
        # Then makes the new config.py
        utils.make_config_file(num_probes, probe_len, probe_dist, min_tm, hairpin_max_tm,
                                offtarget_calc_method, probe_selection_method, filter_repeats)

        print('Altered configuration settings in config.py')
        return()


    def select_file(entry_var, format):
        if format == 'fasta':
            filetypes = (
                ('Fasta file', '*.fa'),
                ('Fasta file', '*.fasta'),
                ('Fasta file', '*.fsa'),
                ('All files', '*.*')
            )
        if format == 'txt':
            filetypes = (
                ('Text file', '*.txt'),
                ('Text file', '*.text'),
                ('Text file', '*.tx'),
                ('All files', '*.*')
            )

        filename = filedialog.askopenfilename(
            title='Open files',
            initialdir='/',
            filetypes=filetypes)
        
        # removes current entry if there is one
        entry_var.delete(0,len(entry_var.get()))
        # puts the file dialog selected file into the entry
        entry_var.insert(0, filename)
        # tries to move the cursor
        # still haven't gotten it to work
        entry_var.xview("end")


    # This initializes the GUI
    root = Tk()
    root.title('COD_FISH Input')

    # create a notebook (the pages)
    notebook = ttk.Notebook(root)
    notebook.pack(pady=10, padx=10, expand=True, fill='both')

    # create and show frames
    main_frame = ttk.Frame(notebook)
    main_frame.pack(fill='both', expand=True)
    advanced_frame = ttk.Frame(notebook)
    advanced_frame.pack(expand=True)

    # add frames to notebook
    notebook.add(main_frame, text='Species and Gene Input')
    notebook.add(advanced_frame, text='Advanced')

    ###################################################
    #       Main frame (Species and Genes)
    ###################################################

    # Species

    # Fill the notebook frame with a LabelFrame
    species_frame = build_label_frame(main_frame, 0, 0, 'Species')

    # Then fill this LabelFrames with input boxes and Labels
    species_label, species_autocombobox = build_input(
        'autocombobox', '', "Please select the species here:", species_frame, 0, 0, 'under', width=35, values=formatted_species)


    # Genes

    gene_frame = build_label_frame(
        main_frame, 0, 3, 'Gene information')

    info_label = ttk.Label(gene_frame, text='You can input gene targets using \nany combination of the following methods below.')
    info_label.grid(column=0, row=3, padx=7, pady=7, sticky='nw')

    # File input
    transcript_file_label, transcript_file_var, gene_symbol_file_button = build_input('filedialog', '', "Transcripts file:", gene_frame, 0, 4, 'under',width=20, format = 'fasta')
    gene_symbol_file_label, gene_symbol_file_var, gene_symbol_file_button = build_input('filedialog', '', "Gene Symbol file:", gene_frame, 0, 7, 'under',width=20, format = 'txt')

    # Symbol and ID input. Can be separated by space, comma, or comma and space
    gene_symbol_label, gene_symbol_var = build_input('entry', '', "Gene symbols:", gene_frame, 0, 10, 'under', width=20)
    gene_ID_label, gene_ID_var = build_input('entry', '', "Ensembl Gene IDs:", gene_frame, 0, 12, 'under', width=20)
    transcript_ID_label, transcript_ID_var = build_input('entry', '', "Ensembl Transcript IDs:", gene_frame, 0, 14, 'under', width=20)

    # Load an image then set it as the window icon.
    # for pyinstaller to make the exe we use this if:

    #icon = PhotoImage(file="codfish_icon.png")
    #root.iconphoto(True, icon)

    # Submit button for running the main part of the program
    button = ttk.Button(main_frame, text='Submit', command=run_program)
    button.grid(column=0, row=13, sticky='se', padx=5, pady=5)
    # makes it so that hitting 'Enter' triggers the submit button.
    button.bind('<Return>', run_program)

    ###################################################
    #       Advanced frame (to set config.py)
    ###################################################

    # Probes

    probes_frame = build_label_frame(advanced_frame, 0, 0, 'Probes')

    probe_num_label, probe_num_entry = build_input(
        'entry', 50, "Number of probes desired per transcript:", probes_frame, 0, 0, 'next_to', width=5, sticky='e')
    probe_len_label, probe_len_entry = build_input(
        'entry', 20, "Length of each probe:", probes_frame, 0, 1, 'next_to', width=5, sticky='e')
    probe_dist_label, probe_dist_entry = build_input(
        'entry', 2, "Minimum distance between each probe:", probes_frame, 0, 2, 'next_to', width=5, sticky='e')
    probe_min_tm_label, probe_min_tm_entry = build_input(
        'entry', 40, "Minimum Tm (\u00b0C) for each probe:", probes_frame, 0, 3, 'next_to', width=5, sticky='e')
    probe_max_hairpin_tm_label, probe_max_hairpin_tm_entry = build_input(
        'entry', 57, "Maximum Tm (\u00b0C) allowed for probe hairpins:", probes_frame, 0, 4, 'next_to', width=5, sticky='e')
    filter_repeats_label, filter_repeats_entry = build_input(
        'combobox', 0, "Filter out probes with nucleotide repeats?:", probes_frame, 0, 5, 'next_to', width=5, sticky='e', values=[True, False])


    # Program parameters

    program_frame = build_label_frame(advanced_frame, 0, 1, 'Program properties')

    offtarget_method_label, off_target_method_entry = build_input(
        'combobox', 0, "Method for calculating \noff target binding:", program_frame, 0, 6, 'next_to', width=20, sticky='e', values=['Tm', 'Alignment_score', 'NUPACK'])
    probe_selection_label, probe_selection_entry = build_input(
        'combobox', 0, "Probe selection method:", program_frame, 0, 7, 'next_to', width=20, sticky='e', values=['Greedy'])

    # Here's the button that triggers the alter_config function above
    advanced_button = ttk.Button(
        advanced_frame, text='Alter Program Parameters', command=alter_config)
    advanced_button.grid(column=0, row=8, sticky='se', padx=5, pady=5)
    # makes it so that hitting 'Enter' triggers the submit button.
    advanced_button.bind('<Return>', alter_config)


    # This tries to see if the program is being run on windows and
    # if so changes a setting which improves the resolution of the GUI
    if os.name == 'nt':
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)

    # This is what constantly updates the window.
    root.mainloop()

# This protects the script from running many times! Not sure why it would do that though
# since it isn't called anywhere else...
if __name__ == "__main__":
    # And this prevents the pyinstaller .exe file from opening the window many times
    # due to incompatibility with multiprocessing
    mp.freeze_support()
    main()
