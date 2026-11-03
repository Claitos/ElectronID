import ROOT                           # pyright: ignore[reportMissingImports]
import uproot                         # pyright: ignore[reportMissingImports]
import XRootD                         # pyright: ignore[reportMissingImports]    
import awkward as ak                  # pyright: ignore[reportMissingImports]
from podio import Frame                         # pyright: ignore[reportMissingImports]
from podio.reading import get_reader            # pyright: ignore[reportMissingImports]

Files=[]

with open('FileList_q2_10to100.txt', 'r') as file:
    lines_list = file.readlines()
    for line in lines_list[:5]: # Read only lines 0:5 - remove [:5] to read all or change 5 to N where N is the number of lines you want
        file_path = line.rstrip() # rstrip to remove trailing white space/new lines   
        Files.append(file_path) # Add file path to array

print("Files to be read:")
for f in Files:
    print(f)
reader = get_reader(Files)


# Use the uproot iterate method to process our list of files - See https://uproot.readthedocs.io/en/stable/uproot.behaviors.TBranch.iterate.html
#for chunk in uproot.iterate({f: "events" for f in Files}, expressions=["MCParticles.PDG"]): # Open files in array f and process events tree with branches specified
    # Process each chunk - Do something
    # print(ak.type(chunk))