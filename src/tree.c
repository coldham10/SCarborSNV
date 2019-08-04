#include <math.h>
#include "tree.h"

int expected_jukes_cantor(long double** JC_dist, long double** freq_numr, int** freq_denom, int m) {
    int i, j;
    long double l_p_bar;
    for (i = 0; i < m; i++) {
        JC_dist[i][i] = logl(0);
        for (j = i + 1; j < m; j++) {
            l_p_bar = freq_numr[i][j] - logl(freq_denom[i][j] * 2);
            l_p_bar = logl(4) - logl(3) + l_p_bar;
            JC_dist[i][j] = JC_dist[j][i] = -(3.0/4.0) * logl(1-expl(l_p_bar));
        }
    }

    return 0;
}

