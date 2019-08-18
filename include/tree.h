/**
 * SCarborSNV: Efficient Phylogeny-aware Single Nucleotide Variant Detection for Single Cells
 *
 * Copyright (C) 2019 Christopher Oldham
 *
 * This file is part of SCarborSNV.
 *
 * SCarborSNV is free software: you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SCarborSNV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SCarborSNV.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef TREE_H_
#define TREE_H_

/*Note in general each node holds information about the edge above it*/
typedef struct Node {
    long double edge_dists[3];
    struct Node* nbrs[3]; /*After rooting: parent[0], left child[1], right child[2]*/
    int id, n_nbrs;
    char  is_root, is_cell;
    /*Upwards step values*/
    long double pi[4];
    long double sum_W_Se, sum_aux0, sum_W3; /*aux0 = p_bar(e)*pi_2(e)/pi_0(e)*/
    /*Downwards step values*/
    long double sum_P_Se;
    long double sum_aux1; /*aux1 = p_bar(e)*pi_1(e)/pi_0(e)*/
    long double partial_prod; /*product of [1-P(Se)] */
    long double sum_W_SL[2]; /*Norm constants passed back up at end of downwards step*/
    /*Node-specific values (not sums, etc)*/
    long double W_Se, aux0, p_bar; /*p_bar is of edge directly above to parent ([0])*/
    long double P_Se, aux1;
    long double W_SL[2]; /*L is current edge, un_normalized partial sum of W(S_e', L_e)*/
    /*Genotype flowdown subproblem objectives*/
    long double P_g[4]; /*Genotype 4 is haploid reference*/
} Node;

    
/* If freq_numr[a][b]/freq_denom[a][b] = p-bar(a,b), JC_dist will
 * be the expected Jukes Cantor distance between pairs of cells */
int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m);
Node* build_tree_nj(long double** dist_mat, int m);
void print_tree(Node* T);
void delete_tree(Node* T);

#endif
