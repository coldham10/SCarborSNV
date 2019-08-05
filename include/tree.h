#ifndef TREE_H_
#define TREE_H_

typedef struct Node {
    long double pi_0, pi_1, pi_2;
    int id;
    struct Node* parent;
    struct Node* l_child;
    struct Node* r_child;
    char  is_root, is_leaf;
} Node;

    
/* If freq_numr[a][b]/freq_denom[a][b] = p-bar(a,b), JC_dist will
 * be the expected Jukes Cantor distance between pairs of cells */
int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m);

#endif
