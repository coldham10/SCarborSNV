import os
import subprocess

def write_fake_bams(m):
    filestxt = open("fakenames.txt", 'w+')
    currdir = os.getcwd()
    for i in range(m):
        fname = currdir + '/fakebams/fake{}.'.format(i)
        samfile = open("fname" + "sam", "w+")
        samfile.write("

        
    print(currdir)

write_fake_bams(0)
