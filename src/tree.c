#include <stdlib.h>
#include <math.h>
#include "tree.h"
/*FIXME remove*/
#include <assert.h>

int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m) {
    int i, j;
    int return_val = 0;
    long double l_p_bar;
    for (i = 0; i < m; i++) {
        JC_dist[i][i] = 0;
        for (j = i + 1; j < m; j++) {
            if (freq_denom[i][j] == 0) {
                return_val = 1;
                JC_dist[i][j] = JC_dist[j][i] = NAN;
                continue;
            }
            l_p_bar = freq_numr[i][j] - logl(freq_denom[i][j] * 2);
            l_p_bar = logl(4) - logl(3) + l_p_bar;
            JC_dist[i][j] = JC_dist[j][i] = -(3.0/4.0) * logl(1-expl(l_p_bar));
        }
    }
    return return_val;
}

void update_r(long double*  r, Node** L, long double** d, int* idx, int n_L);
void min_D(int* out, Node** L, long double** d, long double* r, int* idx, int n_L);
void root_tree(Node* T, Node* parent);

Node* build_tree_nj(long double** dist_mat, int m) {
    /*Following notation in durbin et al. NB assumes full distance matrix*/
    int i, j, last_id;
    int coord[2];
    int n_L = 0;
    Node* node_i;
    Node* root;
    Node** L = malloc(2*m * sizeof(Node*)); 
    /*Indices of active nodes in L*/
    int* L_idx = malloc(2*m * sizeof(int));
    for (i = 0; i < 2*m; i++) { L[i] = NULL; }
    long double* r = malloc(2*m * sizeof(long double));
    /*Extend distance matrix to accommodate internal nodes*/
    long double** d = malloc(2*m * sizeof(long double*));
    for (i = 0; i < 2*m; i++) {
        d[i] = malloc(2*m * sizeof(long double));
        for (j = 0; j < 2*m; j++) {
            d[i][j] = (i <= m && j <= m) ? dist_mat[i][j] : NAN;
        }
    }
    /*Initialize list of avaiable leaves*/
    for (i = 0; i < m + 1; i++) {
        n_L++;
        L_idx[i] = i;
        node_i = malloc(sizeof(Node));
        node_i->id      = i;
        node_i->is_cell = 1;
        node_i->is_root = 0;
        node_i->n_nbrs  = 0;
        node_i->nbrs[0] = node_i->nbrs[1] = node_i->nbrs[2] = NULL;
        L[i] = node_i;
    }
    /*Change root node*/
    L[m]->is_root = 1;
    L[m]->is_cell = 0;
    root = L[m];
    last_id = m;
    /*NJ: Tree building loop */
    while (n_L > 2) {
        update_r(r, L, d, L_idx, n_L);
        min_D(coord, L, d, r, L_idx, n_L);
        /*make new internal node*/
        node_i = malloc(sizeof(Node));
        node_i->id = ++last_id;
        node_i->is_cell = 0;
        node_i->is_root = 0;
        /*Connect neighbours*/
        node_i->nbrs[0] = L[coord[0]];
        node_i->nbrs[1] = L[coord[1]];
        node_i->n_nbrs = 2;
        node_i->nbrs[0]->nbrs[node_i->nbrs[0]->n_nbrs++] = node_i;
        node_i->nbrs[1]->nbrs[node_i->nbrs[1]->n_nbrs++] = node_i;
        /*Compute new distances*/
        node_i->edge_dists[0] = 0.5 * (d[coord[0]][coord[1]] + r[coord[0]] - r[coord[1]]);
        node_i->nbrs[0]->edge_dists[node_i->nbrs[0]->n_nbrs - 1] = node_i->edge_dists[0];
        node_i->edge_dists[1] = d[coord[0]][coord[1]] - node_i->edge_dists[0];
        node_i->nbrs[1]->edge_dists[node_i->nbrs[1]->n_nbrs - 1] = node_i->edge_dists[1];
        /*Remove old nodes from active list*/
        L[coord[0]] = L[coord[1]] = NULL;
        j = 0;
        for (i = 0; i < n_L-2; i++) {
            if (L_idx[i+j] == coord[0] || L_idx[i+j] == coord[1]) {
                j++;
                i--;
                continue;
            }
            L_idx[i] = L_idx[i+j];
        }
        n_L -= 2;
        /*Compute distance of new node to remaining nodes in L*/
        for (i = 0; i < n_L; i++) {
            d[last_id][L_idx[i]] = d[L_idx[i]][last_id] = 0.5 * (d[coord[0]][L_idx[i]] + d[coord[1]][L_idx[i]] - d[coord[0]][coord[0]]);
        }
        /*Add new node to L*/
        L[last_id] = node_i;
        L_idx[n_L++] = last_id;
    }
    /*Join last two*/
    L[L_idx[0]]->nbrs[L[L_idx[0]]->n_nbrs++] = L[L_idx[1]];
    L[L_idx[1]]->nbrs[L[L_idx[1]]->n_nbrs++] = L[L_idx[0]];
    L[L_idx[0]]->edge_dists[L[L_idx[0]]->n_nbrs - 1] = d[L_idx[0]][L_idx[1]];
    L[L_idx[1]]->edge_dists[L[L_idx[1]]->n_nbrs - 1] = d[L_idx[0]][L_idx[1]];
    /*Freeing memory*/
    free(L);
    free(L_idx);
    free(r);
    for (i = 0; i < 2*m; i++) { free(d[i]); }
    free(d);
    root_tree(root, NULL);

    return root;
}

void update_r(long double* r, Node** L, long double** d, int* idx, int n_L) {
    int i, k, a, b;
    for (a = 0; a < n_L; a++) {
        /*Node i is active*/
        i = idx[a];
        r[i] = 0;
        for (b = 0; b < n_L; b++) {
            k = idx[b];
            r[i] += d[i][k];
        }
        /*~average distance between node i and all other active nodes*/
        r[i] /= (n_L - 2);
    }
    return;
}

void min_D(int* out, Node** L, long double** d, long double* r, int* idx, int n_L) {
    int i, j, a, b, min_i, min_j;
    long double D, min_val = INFINITY;
    for (a = 0; a < n_L; a++) {
        /*Node i is active*/
        i = idx[a];
        for (b = a; b < n_L; b++) {
            /*As is j*/
            j = idx[b];
            if (i == j) {
                continue;
            }
            D = d[i][j] - r[i] - r[j];
            if (D < min_val) {
                min_val = D;
                min_i = i;
                min_j = j;
            }
        }
    }
    out[0] = min_i;
    out[1] = min_j;
    return;
}

void root_tree(Node* T, Node* parent) {
    /*Reorganizes data so that the parent is the first neighbor*/
    Node* temp_node;
    long double temp_d;
    int i;
    if (T->is_cell) {
        assert(T->nbrs[0] == parent);
        return;
    }
    for(i = 0; i < 3; i++) {
        if (T->nbrs[i] == parent) {
            break;
        }
    }
    assert(T->nbrs[i] == parent);
    if (i != 0) {
        temp_node = T->nbrs[0];
        T->nbrs[0] = parent;
        T->nbrs[i] = temp_node;
        temp_d = T->edge_dists[0];
        T->edge_dists[0] = T->edge_dists[i];
        T->edge_dists[i] = temp_d;
    }
    root_tree(T->nbrs[1], T);
    if (!T->is_root) {
        root_tree(T->nbrs[2], T);
    }
    return;
}

void delete_tree(Node* T) {
    if (!T->is_cell) {
        delete_tree(T->nbrs[1]);
        if (!T->is_root) {
            delete_tree(T->nbrs[2]);
        }
    }
    free(T);
}



/*
typedef struct Node {
    long double d_0, d_1, d_2;
    long double pi_0, pi_1, pi_2;
    struct Node* nbrs[3]; After rooting: parent, left child, right child
    int id, n_nbrs;
    char  is_root, is_cell;
} Node;
*/
