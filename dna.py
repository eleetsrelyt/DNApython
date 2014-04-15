import random

rna = ["A", "U", "C", "G"]
dna = ["A", "T", "C", "G"]
binary = ["0", "1"]
comp = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    }
tscribe = {
    "A": "U",
    "T": "A",
    "C": "G",
    "G": "C",
    }
tslate = {
    }


def rndout(count=51, array=dna, space=False):
    i = 1
    seq = ""
    if space == True:    
        while i <= count:
            seq += random.choice(array)
            seq += " "
            i += 1
        return seq
    else:
        while i <= count:
            seq += random.choice(array)
            i += 1
        return seq

def complement(seq):
    output = ""
    for i in seq:
        output += comp[i]
    return output

def transcribe(seq):
    output = ""
    for i in seq:
        output += tscribe[i]
    return output
