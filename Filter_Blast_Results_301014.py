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

# Store the Identity by each range
def store_results_2_1(data,temp):
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			# genomic range key
			GK = str(temp[9]) + "-" + str(temp[10])
			
			# Transcript range key
			TK = str(temp[6]) + "-" + str(temp[7])
			
			data[temp[0]][temp[1]][TK] = float(temp[2])
			
		else:
			# genomic range key
			GK = str(temp[9]) + "-" + str(temp[10])
			
			# Transcript range key
			TK = str(temp[6]) + "-" + str(temp[7])
			
			data[temp[0]][temp[1]] = {TK: float(temp[2])}
			
	else:
		# genomic range key
		GK = str(temp[9]) + "-" + str(temp[10])
			
		# Transcript range key
		TK = str(temp[6]) + "-" + str(temp[7])
		
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  {TK: float(temp[2])}
	
	return data


def store_results_3(data,temp):
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			data[temp[0]][temp[1]].append(str(temp[6]) + "-" + str(temp[7]))  
		else:
			data[temp[0]][temp[1]] =  [str(temp[6]) + "-" + str(temp[7])]
	else:
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  [str(temp[6]) + "-" + str(temp[7])]
	
	return data

def store_results_4(data,temp):
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			data[temp[0]][temp[1]].append(str(temp[9]) + "-" + str(temp[10]))  
		else:
			data[temp[0]][temp[1]] =  [str(temp[9]) + "-" + str(temp[10])]
	else:
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  [str(temp[9]) + "-" + str(temp[10])]
	
	return data

# Store Transcript Mapping of Genomic sequences
# e.g. data['CIGGSSA_XXXX']['ssa01']['1233434343-12323232'] = '123-312'
def store_results_5(data,temp):
	if data.has_key(temp[0]):
		if data[temp[0]].has_key(temp[1]):
			# genomic range key
			GK = str(temp[9]) + "-" + str(temp[10])
			
			# Transcript range key
			TK = str(temp[6]) + "-" + str(temp[7])
			
			data[temp[0]][temp[1]][GK] = TK
			
			#data[temp[0]][temp[1]].append(str(temp[9]) + "-" + str(temp[10]))  
		else:
			# genomic range key
			GK = str(temp[9]) + "-" + str(temp[10])
			
			# Transcript range key
			TK = str(temp[6]) + "-" + str(temp[7])
			
			data[temp[0]][temp[1]] = {GK: TK}
			
	else:
		# genomic range key
		GK = str(temp[9]) + "-" + str(temp[10])
			
		# Transcript range key
		TK = str(temp[6]) + "-" + str(temp[7])
		
		data[temp[0]] = {}
		data[temp[0]][temp[1]] =  {GK: TK}
	
	return data

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

def update_coverage_Broken_Genes(BAC_Range,Fos_Range,RID,COV1,IDN,HL,Identity,QC,OF,RANGE,TGM,D_IDN,START_END):
	HITS = []
	M = 0; N = 0
	H = {}
	H1 = {}
	COV_Temp = {}
	for i in BAC_Range:
		#print i
		HITS2 =[]
		temp = i.split("_")
		CID = i.split("_")[0]
		if len(temp) > 2:
			CID = i
		if IDN[RID][CID] >= Identity:
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
			#print RID,i,COVV
			
			# Calculate Identity
			new_IDN = cal_IDN_new(RID,i,RANGE,TGM,D_IDN)
			
			txt = "Single" + "\t" + str(RID) + "\t" + i + "\t" + str(new_IDN) + "\t" + str(COVV*100) + "\n"
			#OF.write(txt)
			# Ex. H["unmapped"] = 0.098
			# Ex. H1["unmapped"] = [`jcf2355567331_0-656_unmapped`]
			
			if COVV != 0:
				if H.has_key(Chr_NO):
					H[Chr_NO] += l2/RL[RID]
					H1[Chr_NO].append(i)
				else:
					H[Chr_NO] = l2/RL[RID]
					H1[Chr_NO] = [i]
						
			M += IDN[RID][CID]	
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
	
	Type = "Combined"	
	
	for i in H1:
		for j in H1[i]:
			CID = j.split("_")[0]
			# Calculate Identity
			new_IDN = cal_IDN_new(RID,j,RANGE,TGM,D_IDN)
			
			if CC*100 >= QC:
				#print Type,"\t",RID,"\t",j,"\t",IDN[RID][j],"\t",COV_Temp[j]*100,"\t",CC*100
				try:
					txt = str(RID) + "\t" + j + "\t" + str(new_IDN) + "\t" + str(COV_Temp[j]*100) + "\t" + str(CC*100) + "\t" + START_END[j] + "\n"
					OF.write(txt)
				except:
					pass
					#txt = str(RID) + "\t" + j + "\t" + str(new_IDN) + "\t" + str(COV_Temp[j]*100) + "\t" + str(CC*100) + "\t" + START_END[j] + "\n"
					#OF.write(txt)
	
	return 0

def cal_IDN_new(Trans,Scaf,RANGE,TGM,D_IDN):
	I = 0
	N = 0
	for i in RANGE[Trans][Scaf]:
		temp = Scaf.split("_")
		N_Scaf = Scaf.split("_")[0]
		if len(temp) > 2:
			N_Scaf = Scaf
		#print temp
		#print Trans,"\t",N_Scaf,"\t",i,"\t",len(temp)
		I += float(D_IDN[Trans][N_Scaf][i])
		
		N += 1
	IDN = I/N
	
	return IDN
	
def fix_overlapping_cov(COV,RANGE,RL,IDN,Identity,QC,HL,OF,D_IDN,TGM,START_END):
	Fos_Range = {}
	Fos_Range_1 = {}
	N = 0
	for i in RANGE:
		N += 1; COVV = 0
		M = 0
		for j in RANGE[i]:
			l,Fos_Range = update_coverage(RANGE[i][j],i,j,Fos_Range,HL)
			CID = j.split("_")[0]
			COV[i][CID] = l/RL[i]
			
			# Calculate Identity
			new_IDN = cal_IDN_new(i,j,RANGE,TGM,D_IDN)
			#print j,"\t",COV[i][CID]*100,"\t",new_IDN
			# Check the coverage and identity 
			if COV[i][CID]*100 >= QC and new_IDN >= Identity:
				CID = j.split("_")[0]
				txt = "ONE\t" + i + "\t" + j  + "\t" + str(new_IDN) + "\t" + str(COV[i][CID]*100) + "\t" + str(COV[i][CID]*100) + "\n"
				#OF.write(txt)
				M = 1
				
			if COV[i][CID] > COVV:
				COVV = COV[i][CID]			
				
		##################################################################
		# Some of the RNA-Seqs are mapping to multiple scaffolds/contigs
		# without overlaps which may mean that these scaffolds/contigs can
		# be joined/put together after further analysis.
		# To list out the list of scaffolds/contigs which may be put together.
		##################################################################
		
		# If the Maximum coverage of a RNA-Seq in any mapping scaffold/Contig
		# is less than 90% then check for mapping scaffolds/contigs of RNA-Seq
		#if COVV*100 < QC:
		M = 0
		if M == 0:
			update_coverage_Broken_Genes(RANGE[i],Fos_Range_1,i,COVV,IDN,HL,Identity,QC,OF,RANGE,TGM,D_IDN,START_END)

	return COV,Fos_Range

def bring_the_match(l,GR):
	for i in GR:
		N = int(i.split("-")[0])
		if l == N:
			match = i
			break
	return match

def sort_genomic_ranges(GR):
	for i in GR:
		for j in GR[i]:
			temp_1 = []
			if len(GR[i][j]) > 1:
				for k in GR[i][j]:
					FN = int(k.split("-")[0])
					temp_1.append(FN)
				temp_1.sort()
				temp_2 = []
				
				# Re-arrange the ranges based on sorted ranges
				for l in temp_1:
					match = bring_the_match(l,GR[i][j])
					temp_2.append(match)
				
				# Replace the OLD ranges by new range
				GR[i][j] = temp_2
				
	return GR

def store_the_new_ranges_to_RANGE(new_ranges,Trans,Scaf,RANGE,M1,TGM):
	Scaf_NEW = Scaf + "_" + str(M1)
	
	NEW_Trans_Ranges = []
	for i in new_ranges:
		NEW_Trans_Ranges.append(TGM[Trans][Scaf][i])
	
	RANGE[Trans][Scaf_NEW] = NEW_Trans_Ranges
		
	return RANGE
def store_start_end(START_END,Range,M1,Scaf):
	if M1 != 0:
		Scaf_NEW = Scaf + "_" + str(M1)
	else:
		Scaf_NEW = Scaf
		
	start = Range[0].split("-")[0]
	end = Range[len(Range)-1].split("-")[1]
	
	START_END[Scaf_NEW] = str(start) + "-" + str(end)
	
	return START_END
	
def pre_process_genomic_range(GR,RANGE,TGM):
	NEW_RANGE = []
	START_END = {}
	for i in GR:
		for j in GR[i]:
			FR = GR[i][j][0]
			End1 = int(FR.split("-")[1])
			N = 0
			temp_range_1 = []
			temp_range_2 = []
			M = 0;M1 = 0
			for k in GR[i][j]:
				if N >= 1:
					End2 = int(k.split("-")[0])
					if (End2 - End1) >= 100000:
						M1 += 1
						START_END = store_start_end(START_END,temp_range_1,M1,j)
						RANGE = store_the_new_ranges_to_RANGE(temp_range_1,i,j,RANGE,M1,TGM)
						temp_range_1 = []
						M = 1
					End1 = int(k.split("-")[1])
					
				temp_range_1.append(k)
				N += 1
				
			if M == 1:
				#print "###:",temp_range_1
				M1 += 1
				START_END = store_start_end(START_END,temp_range_1,M1,j)
				RANGE = store_the_new_ranges_to_RANGE(temp_range_1,i,j,RANGE,M1,TGM)
			else:
				START_END = store_start_end(START_END,temp_range_1,M1,j)
	#for i in START_END:
	#	print i,START_END[i]
	#sys.exit()
	return RANGE,START_END

################################################################
##################### MAIN PROGRAM #############################
################################################################

"""
***********************************************************************
********************** IMPORTANT **************************************
***********************************************************************
New functioans are added to this program due to Torfin's new suggestion.
Removing those functions will give old results

New functions:
sort_genomic_ranges(Genomic_RANGE)
pre_process_genomic_range(Genomic_RANGE,RANGE)

"""

# GO TO THE MAIN DATA	
COV = {}
IDN = {}
D_IDN = {}
RANGE = {}
Genomic_RANGE = {}
Trans_Genomic_Map = {}
RL = {}

# Input BLAST OUTPUT File
In_File = sys.argv[1]

Identity = int(sys.argv[2])

Query_Coverage = int(sys.argv[3])

Hit_Length = int(sys.argv[4])

output_file = sys.argv[5]

OF = open(output_file,"w")

# BLAST OUTPUT FILES
with open(In_File) as infile:
	for i in infile:
		temp = i.split()
		
		# STORE QUERY LENGTH
		RL[temp[0]] = int(temp[8])			
		
		# Store COV
		COV = store_results_1(COV,temp,RL)
		
		# STORE IDN
		IDN = store_results_2(IDN,temp)
		D_IDN = store_results_2_1(D_IDN,temp)
		
		# STORE RANGES
		RANGE = store_results_3(RANGE,temp)

		# STORE RANGES
		Genomic_RANGE = store_results_4(Genomic_RANGE,temp)
		
		# Mapping Transcript ranges to Genomic ranges
		Trans_Genomic_Map = store_results_5(Trans_Genomic_Map,temp)

# Ranges are sorted 
Genomic_RANGE = sort_genomic_ranges(Genomic_RANGE)

RANGE,START_END = pre_process_genomic_range(Genomic_RANGE,RANGE,Trans_Genomic_Map)

##################################################################
# REMOVE OVERLAPPING and Calculate coverage
#################################################################
COV,Fos_Range = fix_overlapping_cov(COV,RANGE,RL,IDN,Identity,Query_Coverage,Hit_Length,OF,D_IDN,Trans_Genomic_Map,START_END)

OF.close()