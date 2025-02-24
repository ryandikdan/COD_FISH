# COD_FISH
![logo_outline](https://user-images.githubusercontent.com/65059714/232628575-97057538-37f3-4d22-b423-b9b5a3f37c53.png)

This repository contains the code for running COD_FISH, the smFISH probe designer. 
## Running the package
To run on windows or mac, we recommend downloading the respective released zip file from https://github.com/ryandikdan/COD_FISH/releases (windows or mac), unzip it, and then double click the gui_COD_FISH.exe (windows) or gui_COD_FISH (mac) file and it should open the GUI and a terminal for progress.
## Running with python
To run it using python, this repository can be downloaded using 
```
git clone https://github.com/ryandikdan/COD_FISH
```
Then we recommend making a virtual environment so that required dependencies don't conflict
```
python -m venv .venv
```
Activate the virtual environment (varies per OS)
```
# Windows Powershell
.\.venv\Scripts\Activate.ps1
# Windows Command Prompt
.venv\Scripts\activate.bat
# Mac or linux
source ./venv/bin/activate
```
You should see an indicator that you're using this environment.
Then install the dependencies in this virtual environment by running
```
pip install -r requirements.txt
```
After this, you should be able to run the GUI via
```
python gui_COD_FISH.py
```
To run via the command line you can use
```
python COD_FISH.py
```
To leave the virtual environment you can simply run
```
deactivate
```

# Next steps

- [ ] Troubleshoot apps and problems across windows and macintosh
- [ ] Make work on M1 and x86 macs
- [ ] Shrink the size of the releases (pyinstaller optimization)
- [ ] Make dynamic programming probe set selection
- [ ] Add more selection criteria (RBP, folding, etc)
- [ ] Rewrite how the program moves around probe sequences and scores (make a class)
- [x] Upload cute logo
- [x] Update to python 3.11
