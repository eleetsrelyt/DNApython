import random

atcg = ["A", "T", "C", "G"]
bin01 = ["0", "1"]

def rndXX(x, array):
    i = 0
    xint = int(x)
    if xint == x and x >= 0:
        while i <= x:
            print random.choice(array),
            i += 1
    else:
        return "You must enter a positive integer integer"        

def rndATCG(x):
    count = 0
    xint = int(x)
    if xint == x and x >= 0:
        while count <= x:
            print random.choice(["A", "T", "C", "G"]),
            count += 1
    else:
        return "You must enter a positive integer integer"

def rnd01(x):
    count = 0
    xint = int(x)
    if xint == x and x >= 0:
        while count <= x:
            print random.choice(["1", "0"]),
            count += 1
    else:
        return "You must enter a positive integer integer"

def rnd101ATCG(x):
    count = 0
    xint = int(x)
    if xint == x and x >= 0:
        while count <= x:
            print random.choice(["1", "0", "A", "T", "C", "G"]),
            count += 1
    else:
        return "You must enter a positive integer integer"
