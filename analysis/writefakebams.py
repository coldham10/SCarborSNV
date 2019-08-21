#!/usr/bin/env python3
import os
import subprocess

def write_fake_bams(m):
    filestxt = open("filenames.txt", 'w+')
    currdir = os.getcwd()
    for i in range(m):
        fname = currdir + '/fakebams/fake{}.'.format(i)
        cell_label = chr(ord('a') + i)
        samfile = open(fname + "sam", "w+")
        samfile.write(("@HD\tVN:1.3\tSO:coordinate\n" +
                "@SQ\tSN:chrM\tLN:16571\n" +
                "@SQ\tSN:chr1\tLN:247249719\n" +
                "@SQ\tSN:chr2\tLN:242951149\n" +
                "@SQ\tSN:chr3\tLN:199501827\n" +
                "@SQ\tSN:chr4\tLN:191273063\n" +
                "@SQ\tSN:chr5\tLN:180857866\n" +
                "@SQ\tSN:chr6\tLN:170899992\n" +
                "@SQ\tSN:chr7\tLN:158821424\n" +
                "@SQ\tSN:chr8\tLN:146274826\n" +
                "@SQ\tSN:chr9\tLN:140273252\n" +
                "@SQ\tSN:chr10\tLN:135374737\n" +
                "@SQ\tSN:chr11\tLN:134452384\n" +
                "@SQ\tSN:chr12\tLN:132349534\n" +
                "@SQ\tSN:chr13\tLN:114142980\n" +
                "@SQ\tSN:chr14\tLN:106368585\n" +
                "@SQ\tSN:chr15\tLN:100338915\n" +
                "@SQ\tSN:chr16\tLN:88827254\n" +
                "@SQ\tSN:chr17\tLN:78774742\n" +
                "@SQ\tSN:chr18\tLN:76117153\n" +
                "@SQ\tSN:chr19\tLN:63811651\n" +
                "@SQ\tSN:chr20\tLN:62435964\n" +
                "@SQ\tSN:chr21\tLN:46944323\n" +
                "@SQ\tSN:chr22\tLN:49691432\n" +
                "@SQ\tSN:chrX\tLN:154913754\n" +
                "@SQ\tSN:chrY\tLN:57772954\n" +
                "@SQ\tSN:chr1_random\tLN:1663265\n" +
                "@SQ\tSN:chr2_random\tLN:185571\n" +
                "@SQ\tSN:chr3_random\tLN:749256\n" +
                "@SQ\tSN:chr4_random\tLN:842648\n" +
                "@SQ\tSN:chr5_random\tLN:143687\n" +
                "@SQ\tSN:chr6_random\tLN:1875562\n" +
                "@SQ\tSN:chr7_random\tLN:549659\n" +
                "@SQ\tSN:chr8_random\tLN:943810\n" +
                "@SQ\tSN:chr9_random\tLN:1146434\n" +
                "@SQ\tSN:chr10_random\tLN:113275\n" +
                "@SQ\tSN:chr11_random\tLN:215294\n" +
                "@SQ\tSN:chr13_random\tLN:186858\n" +
                "@SQ\tSN:chr15_random\tLN:784346\n" +
                "@SQ\tSN:chr16_random\tLN:105485\n" +
                "@SQ\tSN:chr17_random\tLN:2617613\n" +
                "@SQ\tSN:chr18_random\tLN:4262\n" +
                "@SQ\tSN:chr19_random\tLN:301858\n" +
                "@SQ\tSN:chr21_random\tLN:1679693\n" +
                "@SQ\tSN:chr22_random\tLN:257318\n" +
                "@SQ\tSN:chrX_random\tLN:1719168\n" +
                "@RG\tID:{}\tPL:ILLUMINA\tSM:{}\n" +
                "@PG\tID:bwa\tPN:bwa\tVN:0.7.10-r789\n" +
                "@RG\tID:{}\tPL:ILLUMINA\tSM:{}\n").format(*(4*[cell_label])))
        samfile.close()
        subprocess.run(["samtools", "view", "-b", fname + "sam", "-o", fname + "bam"])
        os.remove(fname + "sam")
        filestxt.write(fname + "bam\n")

