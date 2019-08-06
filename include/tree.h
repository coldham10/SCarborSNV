#ifndef TREE_H_
#define TREE_H_

typedef struct Node {
    long double edge_dists[3];
    long double pi_0, pi_1, pi_2;
    struct Node* nbrs[3]; /*After rooting: parent, left child, right child*/
    int id, n_nbrs;
    char  is_root, is_cell;
} Node;

    
/* If freq_numr[a][b]/freq_denom[a][b] = p-bar(a,b), JC_dist will
 * be the expected Jukes Cantor distance between pairs of cells */
int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m);
Node* build_tree_nj(long double** dist_mat, int m);
void print_tree(Node* T);
void delete_tree(Node* T);

#endif
