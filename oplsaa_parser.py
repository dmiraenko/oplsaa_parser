from scipy.sparse import csr_matrix as csr
from argparse import ArgumentParser
from re import search

class Atom():
    symbol = None
    weight = None
    coordiantion = None
    oplsaa_index = None

class OplsaaParser():
    graph = None
    atoms = None

    def _read_hin(self, file : str):
        size = 0
        rstart = [0]
        cols = []
        values = []
        self.atoms = []
        # Get size of graph matrix
        with open(file, "r") as finp:
            while(True):
                while((line := finp.readline()) and not search("^mol \d", line)):
                    pass
                while((line := finp.readline()) and not search("^endmol \d", line)):
                    l = line.strip().split(" ")
                    self.atoms.append(Atom())
                    self.atoms[-1].symbol = l[3]
                    rstart.append(rstart[-1])
                    for n in range(11, 11+2*int(l[10]), 2):
                        cols.append(int(l[n]) - 1)
                        values.append(1)
                        rstart[-1] += 1
                if(not line):
                    break
        self.graph = csr((values, cols, rstart))      

    def __init__(self, file : str):
        if(file.split(".")[-1].lower() == "hin"):
            self._read_hin(file)

    def find_aromatics(self):
        

parser = ArgumentParser(description='Parse molecule files to add OPLSAA coordination data.')
parser.add_argument('-f', '--file', dest='inpfile', default=None, help='location of the input file', metavar="<file>")
args = parser.parse_args()

op_parser = OplsaaParser(args.inpfile)
