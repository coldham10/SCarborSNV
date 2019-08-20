#!/usr/bin/env python3

import numpy as np
import simulator
import csv
import subprocess
import pandas as pd
import os

genotype_dict = {"0/0" : 0,
                "0/1" : 1,
                "1/0" : 1,
                "1/1" : 2,
                "0/." : 0,
                "./0" : 0,
                "1/." : 2,
                "./1" : 2,
                "./." : -1}
def compare_VCF_cells(real_f, inferred_f, n_sites):
    real_f.seek(0)
    inferred_f.seek(0)
    true_reader = list(csv.reader(real_f, delimiter='\t'))
    call_reader = list(csv.reader(inferred_f, delimiter='\t'))
    #Cells
    TP = FN = TN = FP = 0
    for i, field in enumerate(true_reader[0]):
        if field in genotype_dict.keys():
            first_cell_col_t = i
            break
    for i, field in enumerate(call_reader[0]):
        if field in genotype_dict.keys():
            first_cell_col_c = i
            break
    m = len(true_reader[0]) - first_cell_col_t
    for true_line in true_reader:
        found = False
        for call_line in call_reader:
            if call_line[0] == true_line[0] and call_line[1] == true_line[1]:
                #Correct site call found
                found = True
                for i in range(m):
                    if genotype_dict[true_line[first_cell_col_t + i]] > 0:
                        #Cell is really mutant
                        if genotype_dict[call_line[first_cell_col_c + i]] > 0:
                            #And called correctly
                            TP += 1
                        else:
                            #Real mutant not called
                            FN += 1
                    else:
                        #Cell is welltype
                        if genotype_dict[call_line[first_cell_col_c + i]] > 0:
                            #But called as mutant
                            FP += 1
                        else:
                            #Real wt called as wt
                            TN += 1
                break
        if not found:
            # Site mutant but not called
            for i in range(m):
                if genotype_dict[true_line[first_cell_col_t + i]] > 0:
                    #Cell actually mutant
                    FN += 1
                else:
                    TN += 1
    for call_line in call_reader:
        found = False
        for true_line in true_reader:
            if call_line[0] == true_line[0] and call_line[1] == call_line[1]:
                #Already dealt with correctly called sites above
                found = True
                break
        #Site called with no real mutation
        if not found:
            for i in range(m):
                if genotype_dict[call_line[first_cell_col_c + i]] > 0:
                    #Cell called mutant
                    FP += 1
    TN = m * n_sites - TP - FP - FN
    precision = TP/(TP + FP)
    recall    = TP/(TP + FN)
    F1        = (2*TP)/(2*TP + FP + FN)
    return (precision, recall, F1)


def test_phylo(m_cells, iters, params):
    cell_results = pd.DataFrame(columns=["i","Phylo","Precision","Recall","F1"])
    for i in range(iters):
        T = simulator.Phylogeny()
        T.evolve(n_generations=1000, germline=True)
        T.evolve(n_cells=5, germline=False)
        m = len(T.active_nodes)
        T.prepare()
        vcf_f = open("phylo_temp_r.vcf", "w+")
        T.write_vcf(vcf_f)
        vcf_f.close()
        pfile = open("phylo_temp.pileup", "w+")
        T.write_pileup(pfile)
        pfile.close()
        args = ["../SCarborSNV", "-m", str(m),"-p", "phylo_temp.pileup", "-o" "phylo_temp_c.vcf"]
        for name, val in params.items():
            args.append("--{}={}".format(name, val))
        subprocess.run(args)
        real_vcf = open("phylo_temp_r.vcf", "r")
        call_vcf = open("phylo_temp_c.vcf", "r")
        cell_results.loc[2*i] = [i, 'Y', *compare_VCF_cells(real_vcf, call_vcf, 2000)]
        real_vcf.close()
        call_vcf.close()
        args.append("--omit-phlo-inference")
        subprocess.run(args)
        real_vcf = open("phylo_temp_r.vcf", "r")
        call_vcf = open("phylo_temp_c.vcf", "r")
        cell_results.loc[2*i + 1] = [i, 'N', *compare_VCF_cells(real_vcf, call_vcf, 2000)]
        real_vcf.close()
        call_vcf.close()
    os.remove("phylo_temp.pileup")
    os.remove("phylo_temp_r.vcf")
    os.remove("phylo_temp_c.vcf")
    return cell_results

cell_results = test_phylo(10, 30, {})
cell_results.to_csv("phylo_cell_results.csv", index=False)
cell_results = test_phylo(10, 30, {"posterior-threshold": 0})
cell_results.to_csv("phylo_cell_results_no_thresh.csv", index=False)
