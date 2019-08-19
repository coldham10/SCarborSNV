#!/usr/bin/env python3

#  SCarborSNV: Efficient Phylogeny-aware Single Nucleotide Variant Detection for Single Cells
# 
#  Copyright (C) 2019 Christopher Oldham
# 
#  This file is part of SCarborSNV.
# 
#  SCarborSNV is free software: you can redistribute it and/or modify
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  SCarborSNV is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with SCarborSNV.  If not, see <https://www.gnu.org/licenses/>.

import sys

FILES = sys.argv.copy()[1:]
if len(FILES) != 2:
    print("Takes two VCF files: real, called")
    sys.exit()

print(sys.argv)
