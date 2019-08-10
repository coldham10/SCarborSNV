#include <math.h>
#include "inference.h"
#include "math_utils.h"

int upwards_step(Node* T, long double* posteriors) {
    int i;
    long double to_sum[3];
    if (T->is_root) {
        /*Note that information for the root node is stored one below,
         * as each node holds information about the edge directly above it*/
        T->pi[0] = T->pi[1] = T->pi[2] = T->pi[3] = NAN;
        upwards_step(T->nbrs[1], posteriors);
        T->sum_W3 = T->nbrs[1]->sum_W3;
        return 0;
    }
    if (T->is_cell) {
        /*pi_g*/
        for (i = 0; i < 3; i++) {
            T->pi[i] = posteriors[3 * T->id + i];
        }
        /*pi_mu*/
        T->pi[3] = LSE2(T->pi[1], T->pi[2]);
    }
    else {
        /*Make sure cells below are complete*/
        upwards_step(T->nbrs[1], posteriors);
        upwards_step(T->nbrs[2], posteriors);
        /*Update pi_g and pi_mu from children (log space)*/
        for (i = 0; i < 4; i++) {
            T->pi[i] = T->nbrs[1]->pi[i] + T->nbrs[2]->pi[i];
        }
    }
    /*Node specific values*/
    T->p_bar = inverse_JC(T->edge_dists[0]);
    T->W_Se = T->pi[3] + T->p_bar - T->pi[0];
    T->aux0 = T->p_bar + T->pi[2] - T->pi[0];
    /*Upwards partial sums*/
    if (T->is_cell) {
        T->sum_W_Se = T->W_Se;
        T->sum_aux0 = T->aux0;
        T->sum_W3 = T->aux0 + T->p_bar;
    }
    else {
        /*Partial sum of W(S_e)*/
        to_sum[0] = T->nbrs[1]->sum_W_Se;
        to_sum[1] = T->nbrs[2]->sum_W_Se;
        to_sum[2] = T->W_Se;
        T->sum_W_Se = LSE(to_sum, 3);
        /*Partial sum of aux0 (p_bar * pi_2/pi_0)*/
        to_sum[0] = T->nbrs[1]->sum_aux0;
        to_sum[1] = T->nbrs[2]->sum_aux0;
        to_sum[2] = T->aux0;
        T->sum_aux0 = LSE(to_sum, 3);
        /*Overall sum on W^(3)(L)*/
        to_sum[0] = T->nbrs[1]->sum_W3;
        to_sum[1] = T->nbrs[2]->sum_W3;
        to_sum[2] = T->sum_aux0 + T->p_bar;
        T->sum_W3 = LSE(to_sum, 3);
    }
    return 0;
}

int downwards_step(Node* T) {
    static long double total_sum_W_Se;
    int i;
    long double to_sum[3];
    if (T->is_root) {
        total_sum_W_Se = T->nbrs[1]->sum_W_Se;
        T->sum_P_Se = logl(0);
        T->sum_aux1 = logl(0);
        T->partial_prod = logl(1);
        downwards_step(T->nbrs[1]);
        T->sum_W_SL[0] = T->nbrs[1]->sum_W_SL[0];
        T->sum_W_SL[1] = T->nbrs[1]->sum_W_SL[1];
        return 0;
    }
    /*Node specific values*/
    T->P_Se = T->W_Se - total_sum_W_Se;
    T->aux1 = T->pi[1] + T->p_bar - T->pi[0];
    /*Downwards partial sums/product*/ 
    T->sum_P_Se = LSE2(T->P_Se, T->nbrs[0]->sum_P_Se);
    T->sum_aux1 = LSE2(T->aux1, T->nbrs[0]->sum_aux1);
    T->partial_prod = T->nbrs[0]->partial_prod + logl(1-expl(T->P_Se));
    /*More local values*/
    T->W_SL[0] = T->p_bar + T->pi[0] - T->pi[1] + T->sum_aux1;
    T->W_SL[1] = T->p_bar + T->pi[2] - T->pi[1] + T->sum_aux1;
    /*Recurse down*/
    if (T->is_cell) {
        /*No LOH allowed on terminal edges since ADO much more likely*/
        T->sum_W_SL[0] = T->W_SL[0] = logl(0);
        T->sum_W_SL[1] = T->W_SL[1] = logl(0);
        return 0;
    }
    else {
        downwards_step(T->nbrs[1]);
        downwards_step(T->nbrs[2]);
        /*Pass sums back up*/
        for (i = 0; i < 2; i++) {
            to_sum[0] = T->nbrs[1]->sum_W_SL[i];
            to_sum[1] = T->nbrs[2]->sum_W_SL[i];
            to_sum[2] = T->W_SL[i];
            T->sum_W_SL[i] = LSE(to_sum, 3);
        }
        return 0;
    }
}

int DP_genotypes(Node* T, long double* result, long double P_SNV, long double P_LOH) {
    long double P_SNV_e;
    long double P_Le__S[3];
    long double to_sum[4];
    if(T->is_root) {
        /*Root genotype is taken to be definitely 0*/
        T->P_g[0] = logl(1);
        T->P_g[1] = T->P_g[2] = T->P_g[3] = logl(0);
        DP_genotypes(T->nbrs[1], result, P_SNV, P_LOH);
        return 0;
    }
    /*Get global sums from parent*/
    T->sum_W_SL[0] = T->nbrs[0]->sum_W_SL[0];
    T->sum_W_SL[1] = T->nbrs[0]->sum_W_SL[1];
    T->sum_W3      = T->nbrs[0]->sum_W3;
    /*Calculate conditional probabilities*/
    P_SNV_e    = P_SNV + T->P_Se;
    P_Le__S[0] = P_LOH - logl(3) + T->W_SL[0] - T->sum_W_SL[0] - T->sum_P_Se; /*Case 1*/
    P_Le__S[1] = P_LOH - logl(3) + T->W_SL[1] - T->sum_W_SL[1] - T->sum_P_Se; /*Case 2*/
    P_Le__S[2] = P_LOH - logl(3) + T->sum_aux0 + T->p_bar - T->sum_W3; //XXX: should have: - T->partial_prod; /*Case 3*/
    /*Flow down genotype probabilities*/
    /*TODO check these always sum to 1 for each node*/
    /*P(g=0)*/
    to_sum[0] = T->nbrs[0]->P_g[0] + logl(1-expl(P_SNV_e)) + logl(1-expl(P_Le__S[2]));
    to_sum[1] = T->nbrs[0]->P_g[0] + P_SNV_e + P_Le__S[0]; /*SNV and LOH dropped alt on same branch XXX divide by two this term?*/
    to_sum[2] = T->nbrs[0]->P_g[1] + P_Le__S[0];
    T->P_g[0] = LSE(to_sum, 3);
    /*P(g=1)*/
    to_sum[0] = T->nbrs[0]->P_g[0] + P_SNV_e + logl(1-expl(LSE2(P_Le__S[0], P_Le__S[1])));
    to_sum[1] = T->nbrs[0]->P_g[1] + logl(1-expl(LSE2(P_Le__S[0], P_Le__S[1])));
    T->P_g[1] = LSE2(to_sum[0], to_sum[1]);
    /*P(g=2)*/
    to_sum[0] = T->nbrs[0]->P_g[3] + P_SNV_e; /*Silently haploid cell mutates to "homozygous" alt*/
    to_sum[1] = T->nbrs[0]->P_g[0] + P_SNV_e + P_Le__S[1]; /*XXX also divide here by two?*/
    to_sum[2] = T->nbrs[0]->P_g[1] + P_Le__S[1];
    to_sum[3] = T->nbrs[0]->P_g[2];
    T->P_g[2] = LSE(to_sum, 4);
    /*P(g=0 Haploid)*/
    to_sum[0] = T->nbrs[0]->P_g[3] + logl(1-expl(P_SNV_e));
    to_sum[1] = T->nbrs[0]->P_g[0] + P_Le__S[2];
    T->P_g[3] = LSE2(to_sum[0], to_sum[1]);
    if (T->is_cell) {
        result[3*T->id + 0] = LSE2(T->P_g[0], T->P_g[3]);
        result[3*T->id + 1] = T->P_g[1];
        result[3*T->id + 2] = T->P_g[2];
    }
    else {
        DP_genotypes(T->nbrs[1], result, P_SNV, P_LOH);
        DP_genotypes(T->nbrs[2], result, P_SNV, P_LOH);
    }
    return 0;
}




int infer_from_phylogeny(Node* T, long double* posteriors, long double* result, long double P_SNV, long double P_H) {
    upwards_step(T, posteriors);
    downwards_step(T);
    DP_genotypes(T, result, P_SNV, P_H);
    return 0;
}

