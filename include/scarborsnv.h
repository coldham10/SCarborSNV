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
#ifndef SCARBORSNV_H_
#define SCARBORSNV_H_

typedef char nuc_t;

static const nuc_t A = 0;
static const nuc_t C = 1;
static const nuc_t G = 2;
static const nuc_t T = 3;
static const nuc_t INVALID_NUC = 4;
static const nuc_t DELETED_NUC = 5;
static const nuc_t UNKNOWN_NUC = 6; /*For 'N' bases in pileup data */
static const nuc_t PILEEOF_NUC = 7;

typedef struct {
    char mp_fname[512];
    char tmp_fname[512];
    char vcf_fname[512];
    char tree_fname[512];
    double pc_thresh; /*threshold below which cell posterior p(g=0) indicates variant (regardless of phylo)*/
    unsigned int n_threads;
    int m; /*Number of cells */
    int mp_isfile; /*TODO: piling up is to be handled by SCarborSNV */
    int NO_PHYLO; /*Phylogeny aware final probabilities are just posteriors from binomial*/
} global_params_t;

#endif  /*SCARBORSNV_H_ */
