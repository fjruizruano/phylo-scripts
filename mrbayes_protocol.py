#! /usr/bin/python
import os
from subprocess import call
from Bio.Alphabet import IUPAC
from Bio import AlignIO

# data that is introduced by user
file = raw_input("FASTA file name: ")
threads = raw_input("Number of threads: ")
generations = raw_input("Number of generations: ")

try:
	call("mkdir mrbayes", shell=True)
except:
	pass
os.chdir("mrbayes")

# to create nex and nex.tmp files
file_nex = file + ".nex"
file_nex_temp = file + ".nex.temp"

# to convert alignment from FASTA to NEXUS format
ali = AlignIO.read("../" + file, "fasta", alphabet=IUPAC.ambiguous_dna)
AlignIO.write(ali, open(file_nex, "w"), "nexus")
AlignIO.write(ali, open(file_nex_temp, "w"), "nexus")

# to add empty model block to nex.temp for checking alignment
file_nex_data = open(file_nex_temp, "a")
mb_block_test = """BEGIN MRBAYES;
END;
"""

file_nex_data.write(mb_block_test)
file_nex_data.close()

test = True

# to check alignment in MrBayes
try:
	call("mb %s" % file_nex_temp, shell=True)
except:
	test = False

# If nexus file is good, it runs jModelTest
if test == True:
	print "TEST OK! I CONTINUE WITH JMODELTEST.\n"
	os.remove(file_nex_temp) # del nex.temp
	# to run jmodeltest
	call("java -jar /mnt/lmigratoria/jmodeltest-2.1.1/jModelTest.jar -tr %s -d %s -t BIONJ -g 4 -i -f -AIC -w > %s" % (threads, "../"+file, file+".model"), shell=True)
#	print "\nFINISH JMODELTEST\n"

	# to read modeltest results
	model = open(file+".model", "r").readlines()
	model = model[len(model)-1]
	model = model.split("\t")
	model = model[1]

	# dictionaries for parameters of MrBayes for each model
	dic_model = {"GTR":"lset nst=6;\n", "SYM":"lset nst=6;\nprset statefreqpr=fixed(equal);\n", "HKY":"lset nst=2;\n", "K2P":"lset nst=2;\nprset statefreqpr=fixed(equal);\n", "F81": "lset nst=1;\n", "JC":"lset nst=1;\nprset statefreqpr=fixed(equal);\n"}
	dic_gamma_inv = {"I":" rates=propinv", "G":" rates=gamma", "IG":" rates=invgamma"}
	model = model.split("+")
	gaminv = ""
	if len(model) == 3:
		gaminv += dic_gamma_inv["IG"]
	elif len(model) == 2:
		gaminv += dic_gamma_inv[model[1]]
	model_block = dic_model[model[0]]
	model_block_n = model_block.find(";\n")
	model_block = list(model_block) 	
	model_block.insert(model_block_n, gaminv)
	model_block = "".join(model_block)
	print model_block

	# MrBayes block for append to nexus alignment
	mrmodel_block = "BEGIN MRBAYES;\n\t " + "[" + "+".join(model) + "]\n\t" + model_block + """
	mcmcp ngen = %s printfreq=100 samplefreq=100 nchains=4 savebrlens=yes filename=coi;
	mcmc;
	sump filename=coi relburnin=yes;
	sumt filename=coi relburnin=yes conformat=simple;
	END;
	""" % (generations)

	print mrmodel_block

	file_nex_data = open(file_nex, "a")
	file_nex_data.write(mrmodel_block)
	file_nex_data.close()

	# to run created file with MrBayes 
	call("mpirun -np %s mb %s" % (threads, file_nex), shell=True)
