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
#ifndef CALL_VARIANTS_H_
#define CALL_VARIANTS_H_

#include "posteriors.h"

#define CALL_X (-1)
#define CALL_0 (0)
#define CALL_1 (1)
#define CALL_2 (2)

/*Writes a VCF line based on phylogeny aware posteriors and optionally initial posteriors.
 * To only include phylogeny aware calls, set zero_thresh to logl(0)
 * Phylogeny aware calls are made by maximum posterior genotype probability */
int call_to_VCF(FILE* vcf_file, Candidate* c, long double zero_thresh);

#endif
