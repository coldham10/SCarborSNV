#ifndef TREE_H_
#define TREE_H_

typedef struct Node {
    long double d_parent, d_l_child, d_r_child;
    long double pi_0, pi_1, pi_2;
    struct Node* nbr_0; /*After rooting: parent*/
    struct Node* nbr_1; /*After rooting: left child*/
    struct Node* nbr_2; /*After rooting: right child*/
    int id, n_nbrs;
    char  is_root, is_cell;
} Node;

    
/* If freq_numr[a][b]/freq_denom[a][b] = p-bar(a,b), JC_dist will
 * be the expected Jukes Cantor distance between pairs of cells */
int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m);

#endif
