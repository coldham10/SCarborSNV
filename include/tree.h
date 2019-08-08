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
    long double P_silent;
} Node;

    
/* If freq_numr[a][b]/freq_denom[a][b] = p-bar(a,b), JC_dist will
 * be the expected Jukes Cantor distance between pairs of cells */
int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m);
Node* build_tree_nj(long double** dist_mat, int m);
void print_tree(Node* T);
void delete_tree(Node* T);

#endif
