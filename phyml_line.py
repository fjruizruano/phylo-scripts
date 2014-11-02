#!/usr/bin/python

from subprocess import call
from Bio import AlignIO
import os
import sys

try:
	file = sys.argv[1]
except:
	file = raw_input("Introduce FASTA file: ")
try:
	call("mkdir phyml", shell=True)
except:
	pass
os.chdir("phyml")
AlignIO.convert("../"+file, "fasta", file+".phy", "phylip-relaxed")

call("nice phyml -i %s -d nt -b 1000 -b -4 -m GTR -s BEST -v e -c 4 -a e" % (file+".phy"), shell=True)
