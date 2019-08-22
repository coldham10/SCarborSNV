#!/usr/bin/env python3

#  SCarborSNV: Efficient Phylogeny-aware Single Nucleotide Variant Detection for Single Cells
# 
#  Copyright (C) 2019 Christopher Oldham
# 
#  This file is part of SCarborSNV.
# 
#  SCarborSNV is free software: you can redistribute it and/or modify
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  SCarborSNV is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with SCarborSNV.  If not, see <https://www.gnu.org/licenses/>.

import random
import numpy as np

PILEUP_FILENAME = "csimulated_{}.pileup".format("test")
VCF_FILENAME = "csimulated_{}.vcf".format("test")

def to_phred(err):
    Q = np.round(-10 * np.log10(err)) + 33
    if Q > 122:
        Q = 122
    return chr(int(Q))


class Node:

    def __init__(self, num, SNVs={}, LOHs={}, is_root=False):
        self.num = num
        self.SNVs = SNVs
        self.LOHs = LOHs
        self.ADOs = {}
        self.c1 = self.c2 = None
        self.is_root = is_root
        self.print_label = '?'
        return

    def try_mutate(self, tree, germline=False):
        # SNVs
        if random.uniform(0, 1) < tree.SNV_rate:
            mut_site = random.randrange(0, tree.n_sites)
            #Ensure is a new site
            while mut_site in tree.SNVs.keys():
                mut_site = random.randrange(0, tree.n_sites)
            tree.SNVs[mut_site] = random.choice([(0,1), (1,0)])
            if germline and random.uniform(0, 1) < tree.SNV_rate: #HWE
                tree.SNVs[mut_site] = (1,1)
            self.SNVs[mut_site] = tree.SNVs[mut_site]
            return True
        # LOHs
        elif random.uniform(0,1) < tree.haploid_rate:
            LOH_chrom = random.randrange(0, tree.n_chrms)
            if LOH_chrom not in self.LOHs.keys():
                self.LOHs[LOH_chrom] = random.choice([(0,1), (1,0)])
            return True
        else:
            return False


    def try_branch(self, tree):
        if random.uniform(0, 1) < tree.branch_rate:
            self.c1 = Node(tree.last_id+1, SNVs=self.SNVs.copy(), LOHs=self.LOHs.copy())
            self.c2 = Node(tree.last_id+2, SNVs=self.SNVs.copy(), LOHs=self.LOHs.copy())
            return True
        else:
            return False

#N_ITERS = 1000
#CLONAL_ITERS = 1000
#READ_DEPTH = 10
class Phylogeny:
    def __init__(self, n_sites=2000, chrom_length=100, SNV_rate=0.02, branch_rate=0.002, haploid_rate=0.0005, P_ADO_chrom=0.2):
        self.last_id = 0
        self.root = Node(self.last_id, is_root=True);
        self.active_nodes = [self.root]
        self.SNVs = {}
        self.nucs = {}
        self.n_sites        = n_sites
        self.chrom_length   = chrom_length
        self.SNV_rate       = SNV_rate
        self.branch_rate    = branch_rate
        self.haploid_rate   = haploid_rate
        self.P_ADO_chrom    = P_ADO_chrom
        self.n_chrms        = self.n_sites//self.chrom_length
        return

    def evolve(self, n_generations=1000, germline=False, n_cells=0):
        if germline:
            for i in range(n_generations):
                self.root.try_mutate(self, germline=True)
        elif n_cells == 0:
            for i in range(n_generations):
                for cell in self.active_sites.copy():
                    cell.try_mutate(self, germline=False)
                    did_branch = cell.try_branch(self)
                    if did_branch:
                        self.last_id += 2
                        self.active_nodes.remove(cell)
                        self.active_nodes.append(cell.c1)
                        self.active_nodes.append(cell.c2)
        else:
            while len(self.active_nodes) < n_cells:
                for cell in self.active_nodes.copy():
                    cell.try_mutate(self, germline=False)
                    did_branch = cell.try_branch(self)
                    if did_branch:
                        self.last_id += 2
                        self.active_nodes.remove(cell)
                        self.active_nodes.append(cell.c1)
                        self.active_nodes.append(cell.c2)


    def __label_tree(self):
        lab = 'a'
        for cell in self.active_nodes:
            cell.print_label = lab
            lab = chr(ord(lab) + 1)
        return

    def __pick_nucs(self):
        for loc in self.SNVs.keys():
            ref = random.choice(['A','C','G','T'])
            alts = ['A','C','G','T']
            alts.remove(ref)
            alt = random.choice(alts)
            self.nucs[loc] = (ref, alt)
        return

    def __drop_alleles(self):
        for cell in self.active_nodes:
            for c in range(self.n_chrms):
                if np.random.uniform(0,1) < self.P_ADO_chrom:
                    cell.ADOs[c] = random.choice([(0,1), (1,0)])
                    if np.random.uniform(0,1) < self.P_ADO_chrom:
                        cell.ADOs[c] = random.choice([cell.ADOs[c], (1,1)])
        return

    def prepare(self):
        self.__label_tree()
        self.__pick_nucs()
        self.__drop_alleles()


    def __print_node(self, T, f):
        if T.is_root:
            f.write('(')
        if T.c1 == None:
            f.write(T.print_label)
            return
        f.write('(')
        self.__print_node(T.c1, f)
        f.write(":{},".format(1+len(T.c1.SNVs) - len(T.SNVs)))
        self.__print_node(T.c2, f)
        f.write(":{})".format(1+len(T.c2.SNVs) - len(T.SNVs)))
        if T.is_root:
            f.write(":{})RT;\n".format(1 + len(T.SNVs)))
        return

    def print_tree(self, f):
        self.__print_node(self.root, f)
        return


    def write_vcf(self, f):
        for i in range(self.n_sites):
            if i not in self.SNVs.keys():
                continue
            ref, alt = self.nucs[i]
            chrm = i // self.chrom_length
            f.write("Chrm{}\t{}\t.\t{}\t{}\t100\tPASS\t.\tGT".format(
                chrm,
                i % self.chrom_length,
                ref,
                alt))
            for cell in self.active_nodes:
                if i in cell.SNVs.keys():
                    g = list(cell.SNVs[i])
                else:
                    g = [0,0]
                if chrm in cell.LOHs.keys():
                    for pair, val in enumerate(cell.LOHs[chrm]):
                        if val == 1:
                            g[pair] = '.'
                #Uncomment if want VCF to know about ADOs
                #if chrm in cell.ADOs.keys():
                #    for pair, val in enumerate(cell.ADOs[chrm]):
                #        if val == 1:
                #            g[pair] = '.'
                f.write("\t{}/{}".format(*g))
            f.write('\n')
        return


    def __sequence_site(self, locus, init_depth, amp_err=0.002, seq_err=0.02):
        chrm = locus // self.chrom_length
        if locus in self.SNVs.keys():
            ref, alt = self.nucs[locus]
        else:
            ref = random.choice(['A','C','G','T'])
            alts = ['A','C','G','T']
            alts.remove(ref)
            alt = random.choice(alts)
        cells = []
        for cell in self.active_nodes:
            depth = init_depth + int(np.random.normal(0, np.sqrt(init_depth)))
            if depth < 1:
                cells.append((0, '*', '*'))
                continue
            if locus in cell.SNVs.keys():
                g = list(cell.SNVs[locus])
            else:
                g = [0,0]
            if chrm in cell.LOHs.keys():
                for pair, val in enumerate(cell.LOHs[chrm]):
                    if val == 1:
                        g[pair] = '.'
            if chrm in cell.ADOs.keys():
                for pair, val in enumerate(cell.ADOs[chrm]):
                    if val == 1:
                        g[pair] = '.'
            g = [x for x in g if x != '.']
            if len(g) < 1:
                cells.append((0, '*', '*'))
                continue
            reads = []
            quals = []
            for j in range(depth):
                x = ref if np.random.choice(g) == 0 else alt
                if np.random.uniform(0,1) < amp_err:
                    poss = ['A','C','G','T']
                    poss.remove(x)
                    x = np.random.choice(poss)
                curr_seq_err = seq_err + np.random.normal(0, seq_err/5)
                if curr_seq_err < 0.00006:
                    curr_seq_err = 0.00006
                if curr_seq_err > 0.5:
                    curr_seq_err = 0.5
                if np.random.uniform(0,1) < curr_seq_err:
                    poss = ['A','C','G','T','N']
                    poss.remove(x)
                    x = np.random.choice(poss)
                if x == ref:
                    x = '.' if np.random.randint(0,1) == 0 else ','
                reads.append(x)
                quals.append(to_phred(curr_seq_err))
            cells.append((depth, "".join(reads), "".join(quals)))
        return (ref, alt, cells)
            
    def __pileup_write_line(self, f, pos, ref, alt, cells):
        total_depth = 0
        for c in cells:
            total_depth += c[0]
        if total_depth == 0:
            return
        f.write("Chrm{}\t{}\t{}".format(
            pos // self.chrom_length,
            pos %  self.chrom_length,
            ref))
        for cell in cells:
            f.write("\t{}\t{}\t{}".format(
                cell[0],
                cell[1],
                cell[2]))
        f.write("\n")
        return

    def write_pileup(self, f, depth=10, amp_err=0.002, seq_err=0.02):
        for i in range(self.n_sites):
            if i % self.chrom_length == 0:
                current_depth = depth
            elif i % 10 == 0:
                current_depth += int(np.random.normal(0,2))
                if current_depth < 0:
                    current_depth = 0
            ref, alt, cellcalls = self.__sequence_site(i, current_depth, amp_err, seq_err)
            self.__pileup_write_line(f, i, ref, alt, cellcalls)



if __name__ == "__main__":
    T = Phylogeny()
    T.evolve(n_generations=800, germline=True)
    T.evolve(n_cells=10, germline=False)
    print(len(T.active_nodes))
    T.prepare()
    #T.print_tree()
    vcf_f = open(VCF_FILENAME, "w+")
    T.write_vcf(vcf_f)
    vcf_f.close()
    pfile = open(PILEUP_FILENAME, "w+")
    T.write_pileup(pfile)
    pfile.close()

