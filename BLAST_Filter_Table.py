from __future__ import division
from collections import Counter 
import sys
import re
import pickle
import operator

#################################################################
###################### SUB PROGRAMS #############################
#################################################################

# Store the coverage
def store_results_1(data,temp,RL):
	try:
		Cov = float(temp[4])/RL[int(temp[0])]
	except:
		Cov = float(temp[4])/RL[str(temp[0])]
	
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			data[temp[0]][temp[1]] +=  Cov
		else:
			data[temp[0]][temp[1]] =  Cov
	else:
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  Cov
	
	return data

# Store the Identity
def store_results_2(data,temp):
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			data[temp[0]][temp[1]] +=  float(temp[2])
			data[temp[0]][temp[1]] /= 2  
		else:
			data[temp[0]][temp[1]] =  float(temp[2])
	else:
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  float(temp[2])
	
	return data	

def store_results_3(data,temp):
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			data[temp[0]][temp[1]].append(str(temp[7]) + "-" + str(temp[8]))  
		else:
			data[temp[0]][temp[1]] =  [str(temp[7]) + "-" + str(temp[8])]
	else:
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  [str(temp[7]) + "-" + str(temp[8])]
	
	return data

def store_identity(Ind_IDN,temp):
	if int(temp[7]) < int(temp[8]):
		qrange = temp[7] + "-" + temp[8]
	else:
		qrange = temp[8] + "-" + temp[7]
		
	if Ind_IDN.has_key(temp[0]):
		if Ind_IDN[temp[0]].has_key(temp[1]):
			Ind_IDN[temp[0]][temp[1]][qrange] = float(temp[2])
		else:
			Ind_IDN[temp[0]][temp[1]] = {qrange: float(temp[2])}
	else:
		Ind_IDN[temp[0]] = {temp[1]: {qrange: float(temp[2])}}
	
	return Ind_IDN

##########################################################
# COMBINING ALL THE HITS TO CALCULATE THE OVERALL COVERAGE
##########################################################
def update_coverage(BAC_Range,AID,FID,Fos_Range,HT):
	dis = 0
	OL = 0
	HITS = []	
	
	for i in xrange(len(BAC_Range)):		
		temp1 = BAC_Range[i].split("-")
		
		###########################################
		# Length of the hit must be >=100 basepairs
		###########################################
		if (int(temp1[1]) - int(temp1[0])) >= HT or int(temp1[0]) - int(temp1[1]) >= HT:
			if int(temp1[1]) > int(temp1[0]):
				HITS = HITS + range(int(temp1[0]),int(temp1[1]))
			else:
				HITS = HITS + range(int(temp1[1]),int(temp1[0]))
		
	HITS1 = list(set(HITS))
	HITS1.sort()
	l2 = len(HITS1)
	
	
	# Since the hit size should be >=100bps, some hits may not have any results
	if l2 > 0:
		# STORE THE MAPPING RANGE
		if Fos_Range.has_key(AID):
			Fos_Range[AID][FID] = str(HITS1[0]) + "-" + str(HITS1[-1])
		else:
			Fos_Range[AID] = {}
			Fos_Range[AID][FID] = str(HITS1[0]) + "-" + str(HITS1[-1])
	
	return l2,Fos_Range

def update_coverage_Broken_Genes(BAC_Range,Fos_Range,RID,COV1,IDN,HL,Identity,QC,OF):
	HITS = []
	M = 0; N = 0
	H = {}
	H1 = {}
	COV_Temp = {}
	
	for i in BAC_Range:
		HITS2 =[]
		if IDN[RID][i] >= Identity:
			for j in xrange(len(BAC_Range[i])):
				temp1 = BAC_Range[i][j].split("-")
				
				if int(temp1[1]) > int(temp1[0]):
					if int(temp1[1]) - int(temp1[0]) >= HL:
						HITS = HITS + range(int(temp1[0]),int(temp1[1])+1)
				else:
					if int(temp1[0]) - int(temp1[1]) >= HL:
						HITS = HITS + range(int(temp1[1]),int(temp1[0])+1)
					
				if int(temp1[1]) > int(temp1[0]):
					if int(temp1[1]) - int(temp1[0]) >= HL:
						HITS2 = HITS2 + range(int(temp1[0]),int(temp1[1])+1)
				else:
					if int(temp1[0]) - int(temp1[1]) >= HL:
						HITS2 = HITS2 + range(int(temp1[1]),int(temp1[0])+1)
				
			CN = i.split("_")
			Chr_NO = CN[len(CN)-1]
						
			HITS1 = list(set(HITS2))
			HITS1.sort()
			l2 = len(HITS1)
			COVV = l2/RL[RID]
			
			COV_Temp[i] = COVV
			
			# Ex. H["unmapped"] = 0.098
			# Ex. H1["unmapped"] = [`jcf2355567331_0-656_unmapped`]
			
			if COVV != 0:
				if H.has_key(Chr_NO):
					H[Chr_NO] += l2/RL[RID]
					H1[Chr_NO].append(i)
				else:
					H[Chr_NO] = l2/RL[RID]
					H1[Chr_NO] = [i]
						
			M += IDN[RID][i]	
			N += 1
	
	HITS1 = list(set(HITS))
	HITS1.sort()
	l2 = len(HITS1)
	
	CC = l2/RL[RID]
	
	
	##########################################################
	# DECIDE THE TYPE OF THE HITS
	# CASE I: ONLY UNMAPPED, GAPEND, GAPSTART
	# CASE II: ONE CHR WITH UNMAPPED OR GAP START OR GAP END	
	# CASE III: DIFFERENT CHR NUMBERS
	##########################################################
	
	# Check for differnt chromosome numbers
	NDC = 0 # Number of different Chromosomes
	NNC = 0 # Number of NON-chromosomes
	for i in H:
		if re.search("ssa",i):
			NDC += 1
		else:
			NNC += 1
	
	# Check for CASE I
	if NDC == 0 and NNC >= 1:
		Type = "THREE"
	elif NDC == 1 and NNC >= 1:
		Type = "FOUR"
	elif NDC == 1 and NNC == 0:
		Type = "TWO"
	elif NDC >1:
		Type = "FIVE"
	elif NDC == 0 and NNC == 1:
		Type = "TWO"

	################################
	# PRINTING PART
	################################
	Type = "TWO"
	for i in H1:
		for j in H1[i]:
			if CC*100 >= QC:
				#print Type,"\t",RID,"\t",j,"\t",IDN[RID][j],"\t",COV_Temp[j]*100,"\t",CC*100
				#txt = Type + "\t" + str(RID) + "\t" + j + "\t" + str(IDN[RID][j]) + "\t" + str(COV_Temp[j]*100) + "\t" + str(CC*100) + "\n"
				txt =  str(RID) + "\t" + j + "\t" + str(IDN[RID][j]) + "\t" + str(COV_Temp[j]*100) + "\t" + "\n"
				OF.write(txt)				
	
	return 0

# Procduces First outout file for overall plotting
def fix_overlapping_cov_1(COV,RANGE,RL,IDN,Identity,QC,HL,OF1):
	Fos_Range = {}
	Fos_Range_1 = {}
	N = 0
	for i in RANGE:
		N += 1; COVV = 0
		M = 0
		for j in RANGE[i]:
			l,Fos_Range = update_coverage(RANGE[i][j],i,j,Fos_Range,HL)
			COV[i][j] = l/RL[i]
			
			# Check the coverage and identity 
			if COV[i][j]*100 >= QC and IDN[i][j] >= Identity:
				#print "ONE\t", i,"\t",j,"\t",IDN[i][j],"\t",COV[i][j]*100,"\t",COV[i][j]*100
				#txt = "ONE\t" + i + "\t" + j  + "\t" + str(IDN[i][j]) + "\t" + str(COV[i][j]*100) + "\t" + str(COV[i][j]*100) + "\n"
				ID = float("{0:.2f}".format(IDN[i][j]))
				COCC = COV[i][j]*100
				COCC = float("{0:.2f}".format(COCC))
				txt = i + "\t" + j  + "\t" + str(ID) + "\t" + str(COCC) + "\t" + "\n"
				OF1.write(txt)
				M = 1
				
			if COV[i][j] > COVV:
				COVV = COV[i][j]			
				
		##################################################################
		# Some of the RNA-Seqs are mapping to multiple scaffolds/contigs
		# without overlaps which may mean that these scaffolds/contigs can
		# be joined/put together after further analysis.
		# To list out the list of scaffolds/contigs which may be put together.
		##################################################################
		
		# If the Maximum coverage of a RNA-Seq in any mapping scaffold/Contig
		# is less than 90% then check for mapping scaffolds/contigs of RNA-Seq
		#if COVV*100 < QC:
			
		#if M == 0:
		#	update_coverage_Broken_Genes(RANGE[i],Fos_Range_1,i,COVV,IDN,HL,Identity,QC,OF)

	return 0

def choose_idn_group(IDN):
    if IDN >= 98:
        G = ">=98%"
    elif IDN >=95 and IDN < 98:
        G = "95-98%"
    elif IDN >90 and IDN < 95:
        G = "91-94%"
    elif IDN >=81 and IDN <=90:
        G = "81-90%"
    else:
        G = "<=80%"

    return G


# Procduces second outout file for detail plotting
def fix_overlapping_cov_2(COV,RANGE,RL,IDN,Identity,QC,HL,OF2,Ind_IDN):
	Fos_Range = {};Fos_Range_1 = {}
	N = 0
	for i in RANGE:
		COVV = 0;M = 0
		for j in RANGE[i]:
			l,Fos_Range = update_coverage(RANGE[i][j],i,j,Fos_Range,HL)
			COV[i][j] = l/RL[i]
			for k in RANGE[i][j]:
				N += 1
				temp = k.split("-")
				G = choose_idn_group(Ind_IDN[i][j][k])
				if COV[i][j]*100 >= QC and Ind_IDN[i][j][k] >=Identity:
					COVT = float("{0:.2f}".format(COV[i][j]*100))
					txt1 = i + "\t" + j + "\t" + str(Ind_IDN[i][j][k]) + "\t" + str(COVT) + "\t" + temp[0] + "\t" + G + "\t" + str(RL[i]) + "\t" + str(N) + "\n"
					txt2 = i + "\t" + j + "\t" + str(Ind_IDN[i][j][k]) + "\t" + str(COVT) + "\t" + temp[1] + "\t" + G + "\t" + str(RL[i]) +"\t" + str(N) + "\n"
					OF2.write(txt1)
					OF2.write(txt2)
		
	return 0
################################################################
##################### MAIN PROGRAM #############################
################################################################

# GO TO THE MAIN DATA	
COV = {}
IDN = {}
Ind_IDN = {}
RANGE = {}
RL = {}

# Input BLAST OUTPUT File
In_File = sys.argv[1]

Identity = int(sys.argv[2])

Query_Coverage = int(sys.argv[3])

Hit_Length = int(sys.argv[4])

output_file_1 = sys.argv[5]

output_file_2 = sys.argv[6]

OF1 = open(output_file_1,"w")

OF2 = open(output_file_2,"w")

txt = "Query_Seq\tSubject\tIdentity\tQuery_Coverage\n"
OF1.write(txt)

# BLAST OUTPUT FILES
with open(In_File) as infile:
	for i in infile:
		temp = i.split()
		
		# STORE QUERY LENGTH
		RL[temp[0]] = int(temp[13].strip())			
		
		# Store COV
		COV = store_results_1(COV,temp,RL)
		
		# STORE IDN
		IDN = store_results_2(IDN,temp)
		Ind_IDN = store_identity(Ind_IDN,temp)
		
		# STORE RANGES
		RANGE = store_results_3(RANGE,temp)
"""		
for i in Ind_IDN:
	for j in Ind_IDN[i]:
		print i,j,Ind_IDN[i][j]
sys.exit()
"""
##################################################################
# REMOVE OVERLAPPING and Calculate coverage
#################################################################
#COV,Fos_Range = fix_overlapping_cov(COV,RANGE,RL,IDN,Identity,Query_Coverage,Hit_Length,OF)

fix_overlapping_cov_1(COV,RANGE,RL,IDN,Identity,Query_Coverage,Hit_Length,OF1)

fix_overlapping_cov_2(COV,RANGE,RL,IDN,Identity,Query_Coverage,Hit_Length,OF2,Ind_IDN)

OF1.close()
OF2.close()