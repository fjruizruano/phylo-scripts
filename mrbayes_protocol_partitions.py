#! /usr/bin/python
import os
from subprocess import call
from Bio.Alphabet import IUPAC
from Bio import AlignIO

# data that is introduced by user
file = raw_input("FASTA file name: ")
outgroup = raw_input("Taxon as outgroup (empty will not use outgroup): ")
partis = raw_input("Starting positions of each partitions separated by \042-\042: ")
threads = raw_input("Number of threads: ")
generations = raw_input("Number of generations (0 does not run MrBayes): ")

ali = AlignIO.read(open(file), "fasta") 
file_nex = file+".nex" # final file name
file_nex_temp = file+".nex.temp" # temporal file name

# to convert alignment from FASTA to NEXUS format
ali = AlignIO.read(file, "fasta", alphabet=IUPAC.ambiguous_dna)
AlignIO.write(ali, open(file_nex_temp, "w"), "nexus")

# to add empty model block to nex.temp for checking alignment
file_nex_data = open(file_nex_temp, "a")
mb_block_test = """BEGIN MRBAYES;
END;
"""

file_nex_data.write(mb_block_test)
file_nex_data.close()

# to check alignment in MrBayes
try:
	call("mb %s" % file_nex_temp, shell=True)
	print "\nALIGNMENT IS GOOD, NEXT PASS!\n"
	os.remove(file_nex_temp) # del nex.temp
except:
	print "\nPLEASE, CHECK ALIGNMENT FILE.\n\nABORT!\n"
	os.remove(file_nex_temp) # del nex.temp
	sys.exit()

#------------------------------------------------------------------------------

# try to make directory and enter into it
try:
	call("mkdir mrbayes", shell=True)
except:
	pass

os.chdir("mrbayes")

# load partition line
partitions = partis.split("-")

# counter and model memory
i = -1
model_mem = {}

# make partitioned fasta files
for part in partitions:
	i += 1
	start = int(partitions[i])
	try:
		end = int(partitions[i+1])
	except:
		end = ali.get_alignment_length()
	AlignIO.write(ali[:,start:end], open("partition%s.fas" % (str(i)), "w"), "fasta") # split alignment
	call("nice java -jar ~/jmodeltest-2.1.1/jModelTest.jar -tr %s -d %s -t BIONJ -g 4 -i -f -AIC -w > %s" % (threads, "partition"+str(i)+".fas", "partition"+str(i)+".model"), shell=True)
	print "\nFINISH JMODELTEST IN PARTITION\n"

	# to read modeltest results
	model = open("partition"+str(i)+".model", "r").readlines()
	model = model[len(model)-1]
	model = model.split("\t")
	model_mem[i] = model[1]


# -----------------------------------------------------------------------------------

# dictionaries for parameters of MrBayes for each model
dic_model = {"GTR":"    lset nst=6;\n    prset statefreqpr=dirichlet(1,1,1,1);\n", "SYM":"    lset nst=6;\n    prset statefreqpr=fixed(equal);\n", "HKY":"    lset nst=2;\n    prset statefreqpr=dirichlet(1,1,1,1);\n", "K80":"    lset nst=2;\n    prset statefreqpr=fixed(equal);\n", "F81": "    lset nst=1;\n    prset statefreqpr=dirichlet(1,1,1,1);\n", "JC":"    lset nst=1;\n    prset statefreqpr=fixed(equal);\n"}
dic_gamma_inv = {"I":" rates=propinv", "G":" rates=gamma", "IG":" rates=invgamma"}

print model_mem

# begin model block and outgroup
model_block = "BEGIN MRBAYES;\n"

if outgroup != "":
	model_block += "    outgroup %s;\n" % outgroup

# name partitinos (a is added because mrbayes do not reconize numbers)
for a in range(0,len(model_mem)):
	start = int(partitions[a])
	try:
		end = int(partitions[a+1])-1
	except:
		end = ali.get_alignment_length()
	model_block += "    charset " + str(a) + "a\t = %s-%s;\n" % (start, end)

model_block += "    partition favored = %s:" % str(len(model_mem))

for a in range(0,len(model_mem)):
	if a < len(model_mem)-1:
		model_block += " " + str(a) + "a,"
	elif a == len(model_mem)-1:
		model_block += " " + str(a) + "a"

model_block += ";\n    set partition = favored;\n\n"

print model_block

#---------------------------------------------------------------

for a in range(0,len(model_mem)):
	model = model_mem[a]
	model = model.split("+")
	gaminv = ""
	
	# mount model
	if len(model) == 3:
		gaminv += dic_gamma_inv["IG"]
	elif len(model) == 2:
		gaminv += dic_gamma_inv[model[1]]
	model_part = dic_model[model[0]]
	model_part_lset = model_part.find("lset ") + len("lset ")
	model_part_prset = model_part.find("prset ") + len("prset ")
	model_part_salto = model_part.find(";\n")
	model_part = list(model_part)
	model_part.insert (model_part_prset, "applyto=(%s) " % str(a+1))
	model_part.insert(model_part_salto, gaminv)
	model_part.insert (model_part_lset, "applyto=(%s) " % str(a+1))
	model_part.append("\n")
	model_part = "".join(model_part)
	model_block += model_part

print model_block

AlignIO.write(ali, open(file_nex, "w"), "nexus")

# MrBayes block for append to nexus alignment
model_block += """    mcmcp ngen = %s printfreq=100 samplefreq=100 nchains=4 savebrlens=yes filename=coi;
    mcmc;
    sump filename=coi relburnin=yes;
    sumt filename=coi relburnin=yes conformat=simple;
    END;
""" % (generations)
print model_block

file_nex_data = open(file_nex, "a")
file_nex_data.write(model_block)
file_nex_data.close()

# to run created file with MrBayes 
if generations != "0":
	call("nice mpirun -np %s mb %s" % (threads, file_nex), shell=True)
