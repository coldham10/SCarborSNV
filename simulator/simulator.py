#!/usr/bin/env python3

import random
import numpy as np

random.seed(1738)

N_SITES = 10000
CHROM_LENGTH = 1000
HAPLOID_RATE = 0.0001
SNV_RATE = 0.02
BRANCH_RATE = 0.002
N_ITERS = 1400
CLONAL_ITERS = 700
READ_DEPTH = 10

mutation_sites = []
active_nodes   = []
last_id = 0

class Node:

    def __init__(self, num, muts=[], chrms=[], is_root=False):
        self.num = num
        self.mutations = muts
        self.chrmabs = chrms
        active_nodes.append(self)
        self.c1 = self.c2 = None
        self.is_root = is_root
        self.print_label = '?'
        return

    def try_mutate(self):
        if random.uniform(0, 1) < SNV_RATE:
            mut_site = random.randrange(0, N_SITES)
            while mut_site in mutation_sites:
                mut_site = random.randrange(0, N_SITES)
            mutation_sites.append(mut_site)
            self.mutations.append(mut_site)
            return True
        elif random.uniform(0,1) < HAPLOID_RATE:
            n_chrms = N_SITES//CHROM_LENGTH
            #tuple is chrom number, pair number 
            self.chrmabs.append((random.randrange(0, n_chrms), random.randint(0,1)))
        else:
            return False


    def try_branch(self):
        if random.uniform(0, 1) < BRANCH_RATE:
            global last_id
            last_id += 1
            self.c1 = Node(last_id, muts=self.mutations.copy(), chrms=self.chrmabs.copy())
            last_id += 1
            self.c2 = Node(last_id, muts=self.mutations.copy(), chrms=self.chrmabs.copy())
            active_nodes.remove(self)
            return True
        else:
            return False

def label_tree():
    lab = 'a'
    for cell in active_nodes:
        cell.print_label = lab
        lab = chr(ord(lab) + 1)
    return


def print_tree(T):
    if (T.is_root):
        print('(', end='')
    if T.c1 == None:
        print(T.print_label, end='')
        return
    print('(', end='')
    print_tree(T.c1)
    print(":{},".format(len(T.c1.mutations) - len(T.mutations)), end='')
    print_tree(T.c2)
    print(":{})".format(len(T.c2.mutations) - len(T.mutations)), end='')
    if T.is_root:
        print(":{})RT;".format(len(T.mutations)))
    return

def sequence_site():

if __name__ == "__main__":
    root = Node(last_id, is_root=True);
    for i in range(CLONAL_ITERS):
        root.try_mutate()
    for i in range(N_ITERS):
        for n in active_nodes.copy():
            n.try_mutate()
            n.try_branch()
    label_tree()
    print_tree(root)
    current_depth = READ_DEPTH
    for i in range(N_SITES):
        if i % 10 == 0:
            current_depth += np.random.normal(0.1)
