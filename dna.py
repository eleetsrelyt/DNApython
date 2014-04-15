import random

rna = ["A", "U", "C", "G"]
dna = ["A", "T", "C", "G"]
binary = ["0", "1"]
comp = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "N": "N",
    }
tscribe = {
    "A": "U",
    "T": "A",
    "C": "G",
    "G": "C",
    }
tslate3 = {
    'GCA': 'Ala',
    'GCC': 'Ala',
    'GCG': 'Ala',
    'GCT': 'Ala',
    'AGA': 'Arg',
    'AGG': 'Arg',
    'CGA': 'Arg',
    'CGC': 'Arg',
    'CGG': 'Arg',
    'CGT': 'Arg',
    'AAC': 'Asn',
    'AAT': 'Asn',
    'GAC': 'Asp',
    'GAT': 'Asp',
    'TGC': 'Cys',
    'TGT': 'Cys',
    'CAA': 'Gln',
    'CAG': 'Gln',
    'GAA': 'Glu',
    'GAG': 'Glu',
    'GGA': 'Gly',
    'GGC': 'Gly',
    'GGG': 'Gly',
    'GGT': 'Gly',
    'CAC': 'His',
    'CAT': 'His',
    'ATA': 'Ile',
    'ATC': 'Ile',
    'ATT': 'Ile',
    'CTA': 'Leu',
    'CTC': 'Leu',
    'CTG': 'Leu',
    'CTT': 'Leu',
    'TTA': 'Leu',
    'TTG': 'Leu',
    'AAA': 'Lys',
    'AAG': 'Lys',
    'ATG': 'Met',
    'TTC': 'Phe',
    'TTT': 'Phe',
    'CCA': 'Pro',
    'CCC': 'Pro',
    'CCG': 'Pro',
    'CCT': 'Pro',
    'AGC': 'Ser',
    'AGT': 'Ser',
    'TCA': 'Ser',
    'TCC': 'Ser',
    'TCG': 'Ser',
    'TCT': 'Ser',
    'TAA': 'Stop',
    'TAG': 'Stop',
    'TGA': 'Stop',
    'ACA': 'Thr',
    'ACC': 'Thr',
    'ACG': 'Thr',
    'ACT': 'Thr',
    'TGG': 'Trp',
    'TAC': 'Tyr',
    'TAT': 'Tyr',
    'GTA': 'Val',
    'GTC': 'Val',
    'GTG': 'Val',
    'GTT': 'Val',
    }
tslate1 = {
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'AGA': 'R',
    'AGG': 'R',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'AAC': 'N',
    'AAT': 'N',
    'GAC': 'D',
    'GAT': 'D',
    'TGC': 'C',
    'TGT': 'C',
    'CAA': 'Q',
    'CAG': 'Q',
    'GAA': 'E',
    'GAG': 'E',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'CAC': 'H',
    'CAT': 'H',
    'ATA': 'I',
    'ATC': 'I',
    'ATT': 'I',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'TTA': 'L',
    'TTG': 'L',
    'AAA': 'K',
    'AAG': 'K',
    'ATG': 'M',
    'TTC': 'F',
    'TTT': 'F',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'AGC': 'S',
    'AGT': 'S',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TAA': 'Stop',
    'TAG': 'Stop',
    'TGA': 'Stop',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'TGG': 'W',
    'TAC': 'Y',
    'TAT': 'Y',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    }
tslate = {
    'GCA': 'Alanine',
    'GCC': 'Alanine',
    'GCG': 'Alanine',
    'GCT': 'Alanine',
    'AGA': 'Arginine',
    'AGG': 'Arginine',
    'CGA': 'Arginine',
    'CGC': 'Arginine',
    'CGG': 'Arginine',
    'CGT': 'Arginine',
    'AAC': 'Asparagine',
    'AAT': 'Asparagine',
    'GAC': 'Aspartic',
    'GAT': 'Aspartic',
    'TGC': 'Cysteine',
    'TGT': 'Cysteine',
    'CAA': 'Glutamine',
    'CAG': 'Glutamine',
    'GAA': 'Glutamic',
    'GAG': 'Glutamic',
    'GGA': 'Glycine',
    'GGC': 'Glycine',
    'GGG': 'Glycine',
    'GGT': 'Glycine',
    'CAC': 'Histidine',
    'CAT': 'Histidine',
    'ATA': 'Isoleucine',
    'ATC': 'Isoleucine',
    'ATT': 'Isoleucine',
    'CTA': 'Leucine',
    'CTC': 'Leucine',
    'CTG': 'Leucine',
    'CTT': 'Leucine',
    'TTA': 'Leucine',
    'TTG': 'Leucine',
    'AAA': 'Lysine',
    'AAG': 'Lysine',
    'ATG': 'Methionine',
    'TTC': 'Phenylalanine',
    'TTT': 'Phenylalanine',
    'CCA': 'Proline',
    'CCC': 'Proline',
    'CCG': 'Proline',
    'CCT': 'Proline',
    'AGC': 'Serine',
    'AGT': 'Serine',
    'TCA': 'Serine',
    'TCC': 'Serine',
    'TCG': 'Serine',
    'TCT': 'Serine',
    'TAA': 'Stop(Ochre)',
    'TAG': 'Stop(Amber)',
    'TGA': 'Stop(Opal)',
    'ACA': 'Threonine',
    'ACC': 'Threonine',
    'ACG': 'Threonine',
    'ACT': 'Threonine',
    'TGG': 'Tryptophan',
    'TAC': 'Tyrosine',
    'TAT': 'Tyrosine',
    'GTA': 'Valine',
    'GTC': 'Valine',
    'GTG': 'Valine',
    'GTT': 'Valine',
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

def complement(seq=rndout()):
    output = ""
    for i in seq:
        output += comp[i]
    return output

def reverse(seq):
    return seq[::-1]

def transcribe(seq=rndout()):
    output = ""
    for i in seq:
        output += tscribe[i]
    return output

def translate(seq=rndout(), amino=3):
    aaseq = ""
    if amino == 3:
        for i in seq:
            aaseq += tslate3[i]
    return aaseq
