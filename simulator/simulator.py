#!/usr/bin/env python3

import random
import numpy as np

#random.seed(22139)

N_SITES = 2000
CHROM_LENGTH = 100
SNV_RATE = 0.02
BRANCH_RATE = 0.002
HAPLOID_RATE = 0.0005
P_ADO_CHROM = 0.1
AMP_ERR = 0.002
SEQ_ERR = 0.02
N_ITERS = 1000
CLONAL_ITERS = 1000
READ_DEPTH = 10

PILEUP_FILENAME = "csimulated_{}.pileup".format(N_SITES)
VCF_FILENAME = "csimulated_{}.vcf".format(N_SITES)

mutation_sites = []
active_nodes   = []
mut_nucs = {}
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
    print(":{},".format(1+len(T.c1.mutations) - len(T.mutations)), end='')
    print_tree(T.c2)
    print(":{})".format(1+len(T.c2.mutations) - len(T.mutations)), end='')
    if T.is_root:
        print(":{})RT;".format(len(T.mutations)))
    return

def drop_alleles():
    for cell in active_nodes:
        n_chrms = N_SITES//CHROM_LENGTH
        for c in range(n_chrms):
            if np.random.uniform(0,1) < P_ADO_CHROM:
                cell.chrmabs.append((c, np.random.randint(0,1)))
    return

def to_phred(err):
    Q = -10 * np.round(np.log10(err)) + 33
    if Q > 122:
        Q = 122
    return chr(int(Q))

def write_vcf(f):
    for i in range(N_SITES):
        if i not in mutation_sites:
            continue
        ref, alt = mut_nucs[i]
        chrm = i // CHROM_LENGTH
        f.write("Chrm{}\t{}\t.\t{}\t{}\t100\tPASS\t.\tGT".format(
            chrm,
            i % CHROM_LENGTH,
            ref,
            alt))
        for cell in active_nodes:
            if i in cell.mutations:
                if chrm in [x[0] for x in cell.chrmabs]:
                    f.write("\t1/1")
                else:
                    f.write("\t0/1")
            else:
                f.write("\t0/0")
                #TODO let homozygous happen above tree
        f.write('\n')
    return

def pick_nucs():
    for loc in mutation_sites:
        ref = np.random.choice(['A','C','G','T'])
        alts = ['A','C','G','T']
        alts.remove(ref)
        alt = np.random.choice(alts)
        mut_nucs[loc] = (ref, alt)

def sequence_site(locus, init_depth):
    if locus in mutation_sites:
        ref, alt = mut_nucs[locus]
    else:
        ref = np.random.choice(['A','C','G','T'])
        alts = ['A','C','G','T']
        alts.remove(ref)
        alt = np.random.choice(alts)
    cells = []
    for cell in active_nodes:
        depth = init_depth + int(np.random.normal(0, np.sqrt(READ_DEPTH)))
        if depth < 1:
            cells.append((0, '*', '*'))
            continue
        genotype = [0,0]
        if locus in cell.mutations:
            genotype[1] = 1
        chrm = locus//CHROM_LENGTH
        for pos, strand in cell.chrmabs:
            if pos//CHROM_LENGTH == chrm:
                genotype[strand] = -1
        genotype = [g for g in genotype if g != -1]
        if len(genotype) == 0:
            cells.append((0, '*', '*'))
            continue
        reads = []
        quals = []
        for j in range(depth):
            x = ref if np.random.choice(genotype) == 0 else alt
            if np.random.uniform(0,1) < AMP_ERR:
                poss = ['A','C','G','T']
                poss.remove(x)
                x = np.random.choice(poss)
            seq_err = SEQ_ERR + np.random.normal(SEQ_ERR, 0.2)
            if seq_err < 0:
                seq_err = 0
            if np.random.uniform(0,1) < seq_err:
                poss = ['A','C','G','T','N']
                poss.remove(x)
                x = np.random.choice(poss)
            if x == ref:
                x = '.' if np.random.randint(0,1) == 0 else ','
            reads.append(x)
            quals.append(to_phred(seq_err))
        cells.append((depth, "".join(reads), "".join(quals)))
    return (ref, alt, cells)
            
def pileup_write(f, pos, ref, alt, cells):
    total_depth = 0
    for c in cells:
        total_depth += c[0]
    if total_depth == 0:
        return
    f.write("Chrm{}\t{}\t{}".format(
        pos // CHROM_LENGTH,
        pos % CHROM_LENGTH,
        ref))
    for cell in cells:
        f.write("\t{}\t{}\t{}".format(
            cell[0],
            cell[1],
            cell[2]))
    f.write("\n")



if __name__ == "__main__":
    root = Node(last_id, is_root=True);
    for i in range(CLONAL_ITERS):
        root.try_mutate()
    for i in range(N_ITERS):
        for n in active_nodes.copy():
            n.try_mutate()
            n.try_branch()
    label_tree()
    print(len(active_nodes))
    print_tree(root)
    pick_nucs()
    vcf_f = open(VCF_FILENAME, "w+")
    write_vcf(vcf_f)
    vcf_f.close()
    drop_alleles()
    pfile = open(PILEUP_FILENAME, "w+")
    current_depth = READ_DEPTH
    for i in range(N_SITES):
        if i % CHROM_LENGTH == 0:
            current_depth = READ_DEPTH
        elif i % 10 == 0:
            current_depth += int(np.random.normal(0,2)); 
            if current_depth < 0:
                current_depth = 0;
        ref, alt, cellcalls = sequence_site(i, current_depth)
        pileup_write(pfile, i, ref, alt, cellcalls)
    pfile.close()

