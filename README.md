# COD_FISH

![bitmap](https://user-images.githubusercontent.com/65059714/232589678-d63a20d7-d645-4167-9fb0-f9ca5d9bcf40.png)

This repository contains the code for running COD_FISH, the smFISH probe designer. To run on windows or mac, we recommend downloading the respective released zip file from https://github.com/ryandikdan/COD_FISH/releases (windows or mac), unzip it, and then double click the gui_COD_FISH.exe (windows) or gui_COD_FISH (mac) file and it should open the GUI and a terminal for progress. To run it using python, this repository must be downloaded into a directory (and other packages must be installed ie `pip install tqdm primer3-py psutil ensembl_rest`) and running the command `python gui_COD_FISH.py` will open up the GUI or the `python COD_FISH.py` command will run the command line version of the program.

# Next steps

- [ ] Troubleshoot apps and problems across windows and macintosh
- [ ] Shrink the size of the releases (pyinstaller optimization)
- [ ] Make dynamic programming probe set selection
- [ ] Add more selection criteria (RBP, folding, etc)
- [ ] Rewrite how the program moves around probe sequences and scores (make a class)
- [x] Upload cute logo
- [x] Update to python 3.11
