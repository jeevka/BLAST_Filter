from __future__ import division
import sys
import re
import numpy as np
import operator
import pickle

#################################################################
#################### SUB PROGRAMS ###############################
#################################################################
def store_range(RANGE,scaf,qrange,N,HL):
	temp = qrange.split("-")
	
	# Check for the minimum length required 
	if (int(temp[1]) - int(temp[0])) >= HL:
		if RANGE.has_key(scaf):
			RANGE[scaf][N] = qrange  
		else:
			RANGE[scaf] = {1: qrange}
	
	return RANGE

def store_Identity(IDN,scaf,IDN_temp,N):
	if IDN.has_key(scaf):
		IDN[scaf][N] = IDN_temp  
	else:
		IDN[scaf] = {1: IDN_temp}
	
	return IDN


def store_score(score,TScore,scaf):
	if TScore.has_key(scaf):
		TScore[scaf] += float(score)
	else:
		TScore[scaf] = float(score)
	
	return TScore

# This part removes the overlapping hits and
# calculates query coverage
def update_coverage(Range,SCAF,RL,COV,Ind_COV,QID):
	
	HITS = []
	for i in Range:
		temp1 = Range[i].split("-")
		if int(temp1[1]) > int(temp1[0]):
			HITS = HITS + range(int(temp1[0]),int(temp1[1]))
		else:
			HITS = HITS + range(int(temp1[1]),int(temp1[0]))
	
	HITS1 = list(set(HITS))
	HITS1.sort()
	l2 = len(HITS1)
		
        COV[SCAF] = l2/int(RL)
	
	Ind_COV[SCAF] = COV[SCAF] 
	
        return COV,Ind_COV

def fix_overlapping_cov(RANGE,RL,Ind_COV,QID):
    COV = {}
    
    for i in RANGE:
        COV,Ind_COV = update_coverage(RANGE[i],i,RL,COV,Ind_COV,QID)
    
    return COV,Ind_COV

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

def store_start_and_end(QStart,QEnd):
	if int(QStart) < int(QEnd):
		qrange = str(QStart) + "-" + str(QEnd)
	else:
		qrange = str(QEnd) + "-" + str(QStart)
	
	return qrange

# Printing the first part 
def print_print_1(sorted_score,RANGE,Ind_COV,QLength,QID,IDN,Identity,QC,OF1,COV):
	
	if len(sorted_score) >= 5:
	   N = 5
	else: 
	   N = len(sorted_score) 
	
	M1 = 0
	# This Part prints the summary identity and query coverage
	for i in xrange(len(sorted_score)-1,len(sorted_score)-N-1,-1):
	    	ID = sorted_score[i]
	    	M = 1
		T_IDN = 0;N1 = 0
		#print QID,"\t",ID
		if RANGE.has_key(ID):
			for j in RANGE[ID]:
				T_IDN += IDN[ID][M]
				N1 += 1
				M += 1
				M1 += 1
			temp_IDN = T_IDN/N1
			temp_QC = COV[ID]*100
			if temp_IDN >= Identity and temp_QC >= QC:
				#txt1 = QID + "\t" + ID + "\t" + str(T_IDN/N1) + "\t" + str(COV[ID]*100) + "\t" + QLength[QID].strip() + "\n"
				txt1 = QID + "\t" + ID + "\t" + str(T_IDN/N1) + "\t" + str(Ind_COV[ID]*100) + "\t" + QLength[QID].strip() + "\n"
				OF1.write(txt1)

	    	N -= 1
	
	return 0

# Printing the second output file for more detailed line plot
def print_print_2(sorted_score,RANGE,Ind_COV,QLength,QID,IDN,Identity,QC,OF2,LN,COV):
	if len(sorted_score) > 10:
	   N = 10
	else: 
	   N = len(sorted_score)

	M1 = 0;N2 = 0 
	for i in xrange(len(sorted_score)-1,len(sorted_score)-N-1,-1):
	    	ID = sorted_score[i]
	    	M = 1
		T_IDN = 0;#N = 0
		N2 += 1
		if RANGE.has_key(ID):
			for j in RANGE[ID]:
				T_IDN += IDN[ID][M]
				N += 1
				M1 += 1;N1 = 0
				if IDN[ID][M] >= Identity and COV[ID]*100 >= QC:
					LN += 1
					temp = RANGE[ID][j].split("-")
					IDN_G = choose_idn_group(IDN[ID][M])
					txt1 = QID + "\t" + ID + "\t" + str(IDN[ID][M]) + "\t" + str(Ind_COV[ID]*100) + "\t" + str(temp[0]) + "\t" + IDN_G + "\t" + QLength[QID].strip() + "\t" + str(LN) + "\t" + str(N2) +"\n"
					txt2 = QID + "\t" + ID + "\t" + str(IDN[ID][M]) + "\t" + str(Ind_COV[ID]*100) + "\t" + str(temp[1]) + "\t" + IDN_G + "\t" + QLength[QID].strip() + "\t" + str(LN) + "\t" + str(N2) + "\n"
					OF2.write(txt1)
					OF2.write(txt2)
					
				M += 1

	    	N -= 1
	
	return LN

# Separate Identity from the data
def bring_identity(DATA):
	temp = DATA.split()
	t1 = temp[3].replace(")","")
	t2 = t1.replace("(","")
	t3 = t2.replace(",","")
	identity = t3.replace("%","")
	mlen = temp[2].split("/")[1]
	
	return identity


# MAIN SUB-PROGRAM 
def main_main(R1,R2,DATA,G_IDN,QC,HL,OF1,OF2,LN):
	# GO TO THE MAIN DATA
	COV = {}; IDN = {}; RANGE = {}
	Sub_RANGE = {};QLength = {}
	Ind_COV = {}
	
	M = 0;NExpect = 0;FQ = 0; FS = 0
	Score_In = 0;NSeq = 0;NQuery = 0
        NSub = 0;NIdn = 0;TScore = {}
	NScore = 0;NL = 0
	
	# BLAST OUTPUT FILES
	for i in xrange(R1,R2):		
			# Storing the Query Name
			if re.search("Query=",DATA[i]):
				QID = DATA[i].split()[1]
			
			# Storing the length of query
			if re.search("Length=",DATA[i]) and NL == 0:
				RL = DATA[i].split("=")[1]
				NL = 1
				QLength[QID] = RL
			
			# First Sequence 
			if re.search(">",DATA[i]):
				NSeq += 1
				
				# When we reach second sequence we need to store the results again.
				if NSeq > 1:
					# In case if + or - strand
					qrange = store_start_and_end(QStart,QEnd)
					srange = store_start_and_end(SStart,SEnd)

					# Increase the NIdn to store the last result.
					NIdn += 1
					
					# Storing the ranges of the Query match
					RANGE = store_range(RANGE,scaf,qrange,NIdn-1,HL)

					# Storing Subject data also for GBrwose
					Sub_RANGE = store_range(Sub_RANGE,scaf,srange,NIdn-1,HL)
					
					# storing Identities
					IDN = store_Identity(IDN,scaf,IDN_temp,NIdn-1)
				
				temp = DATA[i].split()
				scaf = DATA[i].replace("> ","")				
				scaf = scaf.strip()
				
				NIdn = 0;NQuery = 0;NSub = 0;NScore = 0
				
			# Storing the Score to order the results.
			if re.search("Score =",DATA[i]) and NScore == 0:
				score = DATA[i].split()[2]
				TScore = store_score(score,TScore,scaf)
				NScore = 1

			# Catching Identities
			if re.search("Identities",DATA[i]):
				identity = bring_identity(DATA[i])

				# Store the identity
				IDN_temp = int(identity)
				NIdn += 1
				NQuery = 0
				NSub = 0
				# Store the first identity
				if NIdn == 1:
					IDN = store_Identity(IDN,scaf,IDN_temp,NIdn)
			
			# Catching the first Query start in a hit
			if re.search("Query ",DATA[i]) and NQuery == 0:
				QStart = DATA[i].split()[1]
				NQuery += 1
			
			# Catching every query end to capture the last query end
			if re.search("Query ",DATA[i]):
				QEnd = DATA[i].split()[3]
				
			# When we reach Second Identity within a single sequence
			# We need to store the results.
			if re.search("Identities",DATA[i]) and NIdn > 1:
				# In case if + or - strand
				qrange = store_start_and_end(QStart,QEnd)
				srange = store_start_and_end(SStart,SEnd)

				# Storing the ranges of the Query match
				RANGE = store_range(RANGE,scaf,qrange,NIdn-1,HL)

				# Storing Subject data also for GBrwose
				Sub_RANGE = store_range(Sub_RANGE,scaf,srange,NIdn-1,HL)
				
				# storing Identities
				IDN = store_Identity(IDN,scaf,IDN_temp,NIdn)
	
			# Start and end of of every query and 
			if re.search("Sbjct ",DATA[i]) and NSub == 0:
				SStart = DATA[i].split()[1]
				NSub += 1
				
			if re.search("Sbjct ",DATA[i]):
				SEnd = DATA[i].split()[3]
	
	# In case if + or - strand
	qrange = store_start_and_end(QStart,QEnd)
	srange = store_start_and_end(SStart,SEnd)
	
	# Storing the ranges of the Query match
	RANGE = store_range(RANGE,scaf,qrange,NIdn-1,HL)

        # Storing Subject data also for GBrwose
        Sub_RANGE = store_range(Sub_RANGE,scaf,srange,NIdn-1,HL)				

	# storing Identities
	IDN = store_Identity(IDN,scaf,IDN_temp,NIdn)
		
	# COV and TQC are similar data with similar structure
	COV,Ind_COV = fix_overlapping_cov(RANGE,RL,Ind_COV,QID)
	
	sorted_score = sorted(TScore, key=TScore.get)
	
	# Printing the First part for first file
	print_print_1(sorted_score,RANGE,Ind_COV,QLength,QID,IDN,Identity,QC,OF1,COV)

	#fname = File_name + ".p"
	#pickle.dump(Sub_RANGE, open(fname, "wb" ))
	
	################# This Part prints the Output file for R plot ################	
	LN = print_print_2(sorted_score,RANGE,Ind_COV,QLength,QID,IDN,Identity,QC,OF2,LN,COV)

	return LN

####################################################################################################################################
################################################ MAIN PROGRAM ######################################################################
####################################################################################################################################

# Input File
File_name = sys.argv[1]

# Identity
Identity = int(sys.argv[2])

# Query coverage 
Query_Coverage = int(sys.argv[3])

# Hit length
Hit_Length = int(sys.argv[4])

# Output Files
output_file_1 = sys.argv[5]

output_file_2 = sys.argv[6]

# open the output file 1 
OF1 = open(output_file_1,"w")

txt1 = "Query" + "\t" + "Subject" + "\t" + "Identity" + "\t" + "QueryCov" + "\t" + "QLength\n"
OF1.write(txt1)

# open the output file 2 for R plot
OF2 = open(output_file_2,"w")

# Open the input file
IF = open(File_name,"r")

# Store the ranges of each result in the output file
RESULT_RANGES = []
DATA = IF.readlines()
N = 0
for i in DATA:
	if re.search("Query=",i):
		RESULT_RANGES.append(N)
		
	N += 1
RESULT_RANGES.append(N)

# This part will send the ranges for each result within the file.
LN = 0
for i in xrange(len(RESULT_RANGES)-1):
	# Call the main program
	LN = main_main(RESULT_RANGES[i],RESULT_RANGES[i+1],DATA,Identity,Query_Coverage,Hit_Length,OF1,OF2,LN)
	
# Close the output file
OF1.close()
OF2.close()