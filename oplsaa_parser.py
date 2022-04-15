from multiprocessing.sharedctypes import Value
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
        else:
            raise ValueError("File extension is not supported")

    def find_aromatics(self):
        queue = []
        depth = []
        result = []
        for i in range(len(self.atoms)):
            bondcount = self.graph.indptr[i+1] - self.graph.indptr[i]
            if(self.atoms[i].symbol == "C" and bondcount == 3 or self.atoms[i].symbol == "N" and bondcount == 2):
                queue.append(i)
            try:
                while(queue):
                    if(len(depth) == 6):
                        while(queue[-1] != depth[-1]):
                            if(queue.pop(-1) == depth[0]):
                                result = result + [d for d in depth if d not in result]
                        while(queue[-1] == depth[-1]):
                            queue.pop(-1)
                            depth.pop(-1)

                    while(queue[-1] in depth):
                        queue.pop(-1)
                        while(queue[-1] == depth[-1]):
                            queue.pop(-1)
                            depth.pop(-1)

                    k = queue[-1]
                    depth.append(k)
                    for j in range(self.graph.indptr[k], self.graph.indptr[k+1]):
                        numatm = self.graph.indices[j]
                        if(self.atoms[numatm].symbol == "C" and bondcount == 3 or self.atoms[numatm].symbol == "N" and bondcount == 2):
                            queue.append(numatm)
            except:
                pass
        print(result)

if __name__ == "__main__":
    parser = ArgumentParser(description='Parse molecule files to add OPLSAA coordination data.')
    parser.add_argument('-f', '--file', dest='inpfile', default=None, help='location of the input file', metavar="<file>")
    args = parser.parse_args()

    op_parser = OplsaaParser(args.inpfile)
    op_parser.find_aromatics()
