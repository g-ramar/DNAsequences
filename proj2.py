import random

"""
TO DO LIST:
-[ ] error handle all fasta file input. Make sure it's in actual fasta format
-[ ] method to generate sequence that can function as a riboswitch
-[ ] method to do sequence alignment from protein.fasta to get conserv
-[ ] check if dna sequences are handled and outputed as lowercase or uppercase
-[ ] check if defined functions are utilized (like parse_fasta)
-[ ] add diversity to amino acid sequences
-[ ] change infile to filename for print statements in except IOError statements
-[ ] change lists to immutable lists if we don't need to mutate them??
"""



"""
CONSTANTS are defined here.
"""
AMINO_ACID_TO_CODON = {	'I': ['ATT', 'ATC', 'ATA'],
						'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
						'V': ['GTT', 'GTC', 'GTA', 'GTG'],
						'F': ['TTT', 'TTC'],
						'M': ['ATG'],
						'C': ['TGT', 'TGC'],
						'A': ['GCT', 'GCC', 'GCA', 'GCG'],
						'G': ['GGT', 'GGC', 'GGA', 'GGG'],
						'P': ['CCT', 'CCC', 'CCA', 'CCG'],
						'T': ['ACT', 'ACC', 'ACA', 'ACG'],
						'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
						'Y': ['TAT', 'TAC'],
						'W': ['TGG'],
						'Q': ['CAA', 'CAG'],
						'N': ['AAT', 'AAC'],
						'H': ['CAT', 'CAC'],
						'E': ['GAA', 'GAG'],
						'D': ['GAT', 'GAC'],
						'K': ['AAA', 'AAG'],
						'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
						'STOP': ['TAA', 'TAG', 'TGA']
						}

#promoter/upstream
NEGATIVE_35_SEQUENCE = 'TTGACA' #-35 sequence 'aactgt'
PRIBNOW_BOX = 'TATAAT' #-10 sequence 'ATATTA'
SHINE_DALGARNO = 'AGGAGG'
DEFAULT_EFFECTOR = 'AAUCCUCCGG'
TETRA_LOOP = 'TTCG'
UPSTREAM_LENGTH = []
LENGTHS = []
CHANGABLE= []
LENGTH_OF_X = 0
LOOP_OFFSET = 0
FILLER_SIZE = 0


#terminator
CG_BASEPAIR_MIN = 4 #minimum number of base-pairs for CG rich hairpin loop
POLY_U_TAIL_MIN = 8 #minimum number of sequential Uracil nucleotides

#input filenames
CODON_FREQ_FILENAME = 'codonfreq.txt'
EFFECTOR_FILENAME = 'effector.fasta'
PARAMS_FILENAME = 'params.txt'
PROTEIN_FILENAME = 'protein.fasta'
RESTRIC_ENZ_FILENAME = 'sites.fasta'

#output filename


SITES =  []
"""
CLASSES are defined here
"""
class Seq(object):
	"""basic class to hold a single sequence and its header from a FASTA file.
	This allows access to header/sequence in main method and functions more readable.

	EXAMPLE:
	>>> my_seq = Seq(">This is my header", "ATGCGA")
	>>> my_seq.header
	'>This is my header'
	>>> my_seq.sequence
	'ATGCGA'
	"""
	def __init__(self, header, sequence):
		self.header = header
		self.sequence = sequence

"""
FUNCTIONS are defined here
"""

def print_fasta(seq):
	"""Given a Seq object, print its sequence in FASTA format
	seq - a Seq object
	"""
	assert isinstance(seq, Seq), "Argument seq is not instance of Seq"
	print(seq.header)
	for i in range(0, len(seq.sequence), 80):
		print(seq.sequence[i:i+80])

def write_fasta(seq, outfile):
	"""Given a Seq object, write its sequence in FASTA format to file object
	seq - a Seq object
	outfile - a file object
	"""
	assert isinstance(seq, Seq), "Argument seq is not instance of Seq"
	assert isinstance(outfile, file), "Argument outfile is not instance of file"
	assert outfile.mode == 'w', "Outfile does not have write permission"

	outfile.write("{0}\n".format(seq.header))
	for i in range(0, len(seq.sequence), 80):
		outfile.write("{0}\n".format(seq.sequence[i:i+80]))

def parse_fasta(filename):
	"""Given the file name of a FASTA file, parse the file into a list of Seq objects
	filename - name of file in string

	returns list of Seq objects
	"""

	###TO DO###
	#Error handle incorect fasta files

	infile = None
	try:
		infile = open(filename, 'r')
	except IOError:
		print "Error: Unable to open FASTA file, " + infile
		sys.exit(1)

	#list that will hold Seq objects. this list will be returned at the end of the function
	seq_list = []

	#variable that will temporarily hold the header
	header_holder = ""

	#list that will temporarily hold the sequence. Since FASTA format will only give a max of 80 characters per
	#line, we will store the lines in this list, and then join them together into 1 string. 
	sequence_holder = []

	for line in infile:
		#if line is a new header, then update seq_list with the header_holder and sequence_holder
		if line.startswith(">"):
			#check if header_holder is empty. If this is the first header, then the holder will be empty, so no need to update seq_list.
			if header_holder:
				#update seq_list. Make a Seq object. Use the header stored in header_holder, and use the sequence stored in sequence_holder 
				#(since sequence_holder is a list, join the elements of the list together into 1 string)
				seq_list.append(Seq(header_holder, ''.join(sequence_holder)))
				sequence_holder = []

			#since the current line starts with '>', store this line in header_holder (remove the \n by using rstrip())
			header_holder = line.rstrip()

		#keep adding to sequence_holder if line is not a header
		else:
			#since the line is not a header, it's a sequence. Append this line to sequence_holder after removing \n via rstrip()
			sequence_holder.append(line.rstrip())

	#for loop doesn't update the last protein sequence, so update it manually
	seq_list.append(Seq(header_holder, ''.join(sequence_holder)))

	#error handle here

	return seq_list


def parse_params(filename):
	"""Given the file name of a parameters file, parse the file for parameter N.

	filename - file name of parameters file

	returns parameter N
	"""

	"""
	TO DO:
	-use parse fasta
	"""
	infile = None
	try:
		infile = open(filename, 'r')
	except IOError:
		print ("Error: Unable to open Parameters file, " + filename)
		sys.exit(1)

	#pulls the lines
	count = 0
	for line in infile:
                count += 1
		if count == 1:
                    try:
                        value = int(line)
                    except ValueError:
                        print("Error: Parameters not castable to int, " + line)
                        sys.exit(1)
		elif (count >= 1 and line != None):
		    raise IOError("Error: Parameter file not properly formed, " + line)
                    sys.exit(1)
	return value
    
N_variant = parse_params(PARAMS_FILENAME)


def parse_sites(filename):
	"""Given the file name of restriction enzyme sites, parse the file into a list of sequences.

	filename - file name of restriction enzyme sites

	returns list of sequences
	"""

	"""
	TO DO:
	-use parse fasta
	"""
	infile = None
	try:
		infile = open(filename, 'r')
	except IOError:
		print "Error: Unable to open Restriction Sites file, " + infile
		sys.exit(1)

	header = ''
	sites = []
	#adds all the sites to a particular sites array
	for line in infile:
		if line.startswith(">"):
			header = line.rstrip()
		else:
			sites.append(line.rstrip())
	
	#computes the reverse complement for each site and also adds it into the sites array
	temp_sites = []
	for site in sites:
		temp_sites.append(reverse_complement(site))

	#adds temp_sites to sites
	sites += temp_sites

	return sites

def reverse_complement(sequence):
	"""Given a DNA sequence, generate its reverse complement.

	sequence - DNA sequence (string)

	returns reverse complement of sequence (string)
	"""

	reverse_seq = sequence[::-1]
	complement = ''
	for x in reverse_seq:
		if x == "A":
			complement += 'T'
		elif x == 'T':
			complement += 'A'
		elif x == "C":
			complement += 'G'
		elif x == "G":
			complement += 'C'
	return complement



def validseq(sequence):
	for site in SITES:
		if dna_sequence.find(site) != -1:
			return False
	return True


def parse_codonfreq(filename):
	"""Given the file name of a codon frequency txt file, parse the file into a dictionary
	with key-value pairs of codons to their new weighted codon frequency.

	filename - name of file (string)

	returns weighted codon frequencies (dictionary)
	"""

	###TO-DO###
	#error handle

	infile = None
	try:
		infile = open(filename, 'r')
	except IOError:
		print "Error: Unable to open Codon Frequency file, " + infile
		sys.exit(1)

	raw_codon_freqs = {} #dictionary with all codons and its frequency {codon:frequency}
	weighted_codon_freqs ={}#dictionary with codons and weighted frequency relative to the aa they represent {codon: new frequency}

	for line in infile:
		line = line.split() #turns 'ATC 0.25' into ['ATC', '0.25']
		codon_name = line[0] #index 0 is the codon sequence
		codon_freq = float(line[1]) #index 1 is the frequency, but we want the numerical value, not the string

		raw_codon_freqs[codon_name] = codon_freq #add this codon and its frequency to dictionary

	for aa in AMINO_ACID_TO_CODON:
		total = 0.0
		aa_codon_list = AMINO_ACID_TO_CODON[aa] #gets all codons (in a list) that represent current amino acid
		for codon in aa_codon_list:
			#add this codon's freq to total
			total += raw_codon_freqs[codon.lower()]

		#go through all codons again, this time assign the new codon freq to weighted_codon_freqs
		for codon in aa_codon_list:
			weighted_codon_freqs[codon.lower()] = raw_codon_freqs[codon.lower()]/total

	return weighted_codon_freqs

def generate_aa_seq(seq_list):
	"""Given an alignment, generates an amino acid sequence by randomly choosing an amino acid from the alignment at every index position.

	seq_list - a list of seq objects
	
	returns amino acid sequence (string)
	"""

	def homologybuilderfinalmega(objectlist):
		listnonpolar = ["A", "G", "I", "L", "V", "F", "W"]
		listpolar = ["Y", "N", "C", "Q", "M", "S", "T", "P"]
		acidic = ["D", "E"]
		basic = ["R", "H", "K"]
		allacids = listnonpolar + listpolar + acidic + basic

		def homologybuilderfinal(listobjects):
			listsequences = [x.sequence for x in listobjects]
			homologylist = []
			for i in range(len(listsequences)):
				for j in range(len(listsequences[i])):
					try:
						homologylist[j] += listsequences[i][j]
					except:
						homologylist.append(listsequences[i][j])
			return homologylist

		def counterpercent(string, array):
			#print("counterpercent string: " + string)
			#print("counterpercent array: " + str(array))
			count = 0.0
			for x in array:
				count += string.count(x)
			#print("counterpercent: " + str(count))
			#print("counterpercent return: " + str(count/len(string)))
			return count/len(string)

		def countertotal(string, value = allacids):
			maxvalue = 0.0
			for x in value:
				maxvalue = max(maxvalue, string.count(x))
			return maxvalue

		def countall(string, value = allacids):
			count = 0
			for x in value:
				if x in string:
					count += 1
			return count

		def homologybuilderfinal2(homologylist):
			homologyliststats = []
			for x in homologylist:
				homologyliststats.append([x, [countall(x),
										  countertotal(x) / len(x),
										  counterpercent(x, basic),
										  counterpercent(x, acidic),
										  counterpercent(x, listnonpolar),
										  counterpercent(x, listpolar)]])
			return homologyliststats

		"""homologyliststats of form [['string', [statsvalues]], ...]"""
		def homologybuilderfinal3(homologyliststats):
			finalsequence = ""
			for x in homologyliststats:
				if x[1][0] < 6 or x[1][1] > .50:
					finalsequence += random.choice(x[0])
				else:
					basicval = x[1][2]
					acidicval = basicval + x[1][3]
					nonpolarval = acidicval + x[1][4]
					polarval = nonpolarval + x[1][5]
					randomnum = random.random()

					#print(x)

					if randomnum <= basicval:
						#print("basicval")
						finalsequence += random.choice(basic)
					elif randomnum <= acidicval:
						#print("acidicval")
						finalsequence += random.choice(acidic)
					elif randomnum <= nonpolarval:
						#print("nonpolarval")
						finalsequence += random.choice(listnonpolar)
					elif randomnum <= polarval:
						#print("else")
						finalsequence += random.choice(listpolar)
			return finalsequence

		return homologybuilderfinal3(homologybuilderfinal2(homologybuilderfinal(objectlist)))

	return homologybuilderfinalmega(seq_list)


def convert_aa_to_dna(aa_sequence, codon_freq):
	"""Given an amino acid sequence and the amino acid codon frequencies, 
	generates the dna sequence that codes for the aa sequence.

	aa_sequence - amino acid sequence (string)
	codon_freq - codons and their (weighted) frequencies (string)

	returns dna sequence (string)
	"""
	proteinsequence = aa_sequence

	def convert_aa_to_codon(aa):
		"""converts a given amino acid to a codon by using the codon frequency table
		aa - amino acid (string)

		returns codon (string)
		"""
		listof_codons = AMINO_ACID_TO_CODON[aa]

		upto = 0
		for codon in listof_codons:
			rand_num = random.random()
			upto = upto + codon_freq[codon.lower()]
			if rand_num < upto:
				return codon

		#error handling, code shoudn't reach here since a codon should be returned in the for loop
		assert False, "convert_aa_to_codon failed, codon not returned"

	def shuffler(x, n = 1000):
		count = 0
		def listwrong(x2):
			for i in range(len(x2) - 2):
				if x2[i] == x2[i+1]:
					return True
			return False

		while (listwrong(x) and count < n):
			count += 1
			value = x.pop(0)
			index = int(random.random() * len(x))
			x.append(x[index])
			x[index] = value
		return x
		
	G = shuffler([convert_aa_to_codon("G") for x in range(proteinsequence.count('G'))])
	A = shuffler([convert_aa_to_codon("A") for x in range(proteinsequence.count('A'))])
	V = shuffler([convert_aa_to_codon("V") for x in range(proteinsequence.count('V'))])
	L = shuffler([convert_aa_to_codon("L") for x in range(proteinsequence.count('L'))])
	I = shuffler([convert_aa_to_codon("I") for x in range(proteinsequence.count('I'))])
	M = shuffler([convert_aa_to_codon("M") for x in range(proteinsequence.count('M'))])
	F = shuffler([convert_aa_to_codon("F") for x in range(proteinsequence.count('F'))])
	W = shuffler([convert_aa_to_codon("W") for x in range(proteinsequence.count('W'))])
	P = shuffler([convert_aa_to_codon("P") for x in range(proteinsequence.count('P'))])
	S = shuffler([convert_aa_to_codon("S") for x in range(proteinsequence.count('S'))])
	T = shuffler([convert_aa_to_codon("T") for x in range(proteinsequence.count('T'))])
	C = shuffler([convert_aa_to_codon("C") for x in range(proteinsequence.count('C'))])
	Y = shuffler([convert_aa_to_codon("Y") for x in range(proteinsequence.count('Y'))])
	N = shuffler([convert_aa_to_codon("N") for x in range(proteinsequence.count('N'))])
	Q = shuffler([convert_aa_to_codon("Q") for x in range(proteinsequence.count('Q'))])
	D = shuffler([convert_aa_to_codon("D") for x in range(proteinsequence.count('D'))])
	E = shuffler([convert_aa_to_codon("E") for x in range(proteinsequence.count('E'))])
	K = shuffler([convert_aa_to_codon("K") for x in range(proteinsequence.count('K'))])
	R = shuffler([convert_aa_to_codon("R") for x in range(proteinsequence.count('R'))])
	H = shuffler([convert_aa_to_codon("H") for x in range(proteinsequence.count('H'))])
	upstreamsequence = ""
	for x in proteinsequence:
		value = eval(x).pop(0)
		upstreamsequence += value
	return upstreamsequence

def makefixedsites(output_seq):

	"""
	FORMAT:										str index 			length
	-35 seq (index -35 to -30) 					(0 to 5)			(6)
	filler (index -29 to -11)					(6 to 24)			(19)
	-10 seq (index -10 to -5) 					(25 to 30)			(6)
	filler (index -4 to -1)						(31 to 34)			(4)
	START TRANSCRIPTION							(35 to X) 			(x - 35 + 1)
	shine dalgarno seq 							(X + 1 to X + 6)	(6)
	filler										(X + 7 to X + 8)	(2)
	START CODON CODING REGION					(X + 9 ....)		(...)
	"""
	#fill in here to implement a way when combining sequences we keep track of the offset & store which sites we need to keep fixed
	#following example output, start of transcription site = 0, 40 basepairs before it, so it will now have an offset of 40

	#iterate through the indices of the output sequence
	for i in range(len(output_seq)):
		#sequence is -35 sequence that cannot change
		if (i >= 0 and i <= 5):
			CHANGEABLE[i] = 1
		#sequence is -10 sequence that cannot change
		elif (i >= 25 and i <=30):
			CHANGEABLE[i] = 1
		#the riboswitch that follows the start of the transcription site
		elif (i >= 35 and i <= 44 + LENGTH_OF_X + LENGTHS[1]):
			#for bases in the loop, these can vary since they do not base, and therefore we set CHANGEABLE of these indices to 0
			if(i >= 35  and i <= 35 + LENGTH_OF_X):
				if(i >= 35 + LOOP_OFFSET  and i < 35 + LOOP_OFFSET + LOOP_LENGTH):
					CHANGEABLE[i] = 0
				#everything else in the riboswitch cannot be changed for the effector to bind properly
				else:
					CHANGEABLE[i] = 1
			#Shine Dalgarno sequence that must br present in the sequence
			elif(i> 35+LENGTH_OF_X and i <= 41+ LENGTH_OF_X):
				CHANGEABLE[i] = 1
			#filler bases that we can arbitrarily change if there is a restriction enzyme/reverse complement present in the output_seq
			elif(i > 41+LENGTH_OF_X and i <= 43 + LENGTH_OF_X):
				CHANGEABLE[i] = 0
			#start of the DNA coding region
			else:	
				#if the index is not the start of a codon then we cannot change it. (i - 44 - LENGTH_OF_X) is the offset of the CHANGEABLE array from the start of the DNA_SEQ
				if((i-44-LENGTH_OF_X) %3 != 0):
					CHANGEABLE[i] = 1
				else:
					CHANGEABLE[i] = 0
				
		else:
			CHANGEABLE[i] = 0
			#filler sequence before the terminator sequence that we can arbitrarily change by basepair
			#if(i <= 44+LENGTH_OF_X+LENGTHS[1] + FILLER_SIZE):
				#CHANGEABLE[i] = 0
			#terminator sequence that we cannot change
			#else:
				#CHANGEABLE[i] = 1

def checkconstraints(output_seq,aa_sequence):
	#Store the array of all of the fixed sites of the output sequence
	SITES_NEEDED = makefixedsites(output_seq)

	#Keep changing the sequence until we find a valid sequence
	while(validseq(output_seq) == False):
		#for all of the restriction sites, check where the restriction strict occurs, and store that in the variable index
		for site in SITES:
			index = output_seq.index(site)
			#if the CHANGEABLE list suggests that we cannot change the certain base, lets keep incrementing the index
			while(CHANGEABLE[index] == 1):
				index +=  1
			#if the index is within the coding region of the DNA, we have to compute the offset from the index to the start of the DNA sequence
			if(index > 43 + LENGTH_OF_X and i <= 43 + LENGTH_OF_X + LENGTHS[1]):
				index_offset = index - 44 - LENGTH_OF_X
				aa_seq_index = index_offset / 3
				aa_to_replace = aa_sequence[aa_seq_index]
				#the codon that we want to replace generated by a call to convert_aa_to_codon for the particular amino acid
				new_codon = convert_aa_to_codon(aa_to_replace)
				#reassigning the bases of the output_sequence to reflect the new codon
				output_seq[index] = new_codon[0]
				output_seq[index+1] = new_codon[1]
				output_seq[index+2] = new_codon[2]
			#index not in the coding region, we can replace it base by base, so from the list of possible bases, eliminate the current base, and randomly pick a new base
			else:
				bases_list = ['A', 'C','T','G']
				bases_list = bases_list.remove(output_seq[index])
				output_seq[index] = bases_list[random.randint(0,2)]








def makevalidDNAseq(dna_sequence,aa_sequence,aa_codon_freq):
	while validDNAseq():
		for site in SITES:
			index = dna_sequence.index(site)
			while(index % 3 != 0):
				index += 1
			aa_seq_index = index / 3;
			aa_to_replace = aa_sequence[aa_seq_index]
			new_codon = convert_aa_to_codon(aa_to_replace,aa_codon_freq)
			dna_sequence[index] = new_codon[0]
			dna_sequence[index+1] = new_codon[1]
			dna_sequence[index+2] = new_codon[2]

	return dna_sequence

def DNAstrandtoRNA(dna_sequence):
	rna_sequence = ''
	for letters in dna_sequence:
		if letters == 'T':
			rna_sequence += 'U'
		else:
			rna_sequence += letters

	return rna_sequence

def generate_upstream(filename):
	"""Given the file name for the effector, generates the promoter sequence and a functional riboswitch sequence.

	filename - file name of effector file (string)

	returns upstream sequence (string)
	"""

	"""
	FORMAT:										str index 			length
	-35 seq (index -35 to -30) 					(0 to 5)			(6)
	filler (index -29 to -11)					(6 to 24)			(19)
	-10 seq (index -10 to -5) 					(25 to 30)			(6)
	filler (index -4 to -1)						(31 to 34)			(4)
	START TRANSCRIPTION							(35 to X) 			(x - 35 + 1)
	shine dalgarno seq 							(X + 1 to X + 6)	(6)
	filler										(X + 7 to X + 8)	(2)
	START CODON CODING REGION					(X + 9 ....)		(...)
	"""

	"""TO DO:
	UPPERCASE VS LOWERCASE CHECK
	EFFECTOR IS IN RNA NOT DNA
	"""
	def generate_riboswitch(effector):
		index_a = effector.find('a')
		temp = effector[::-1].find('a')

		if temp < index_a:
			index_a = temp
			effector = effector[::-1]

		effector_in_stem = effector[0:index_a]
		remaining_effector = effector[index_a:]

		stem_strand = generate_random_seq(random.randint(10,20), True) + effector_in_stem
		revcomp_stem_strand = reverse_complement(stem_strand)
		loop = TETRA_LOOP #or generate_random_loop(random.int(5, 8))

		revcomp_remaining_effector = reverse_complement(remaining_effector)

		riboswitch_sequence = revcomp_remaining_effector + revcomp_stem_strand + loop + stem_strand

		return riboswitch_sequence

	def generate_alternative_riboswitch(effector):
		edge_base = SHINE_DALGARNO[0]
		comp_base = ''
		if edge_base == 'A':
			comp_base = 'T'
		elif edge_base == 'T':
			comp_base = 'A'
		elif edge_base == 'C':
			comp_base = 'G'
		elif edge_base == 'G':
			comp_base = 'C'
		else:
			assert False, 'base-pair error'

		"""RNA TO DNA"""

		revcomp_effector = reverse_complement(effector)

		filler = generate_random_seq(random.randint(10, 15))
		revcomp_filler = reverse_complement(filler)

		loop = generate_random_loop(random.randint(4,8))

		index_to_remove = [0, len(effector)/2, len(effector)/4]

		bases = ['A', 'T', 'C', 'G']

		dna_effector = effector

		for i in index_to_remove:
			base = dna_effector[i]
			rand_base = bases.choice()
			while rand_base == base:
				rand_base = bases.choice()

			dna_effector = dna_effector[0:i] + rand_base + dna_effector[i + 1:]

		first_filler = ''

		for i in range(1, len(SHINE_DALGARNO)):
			rand_base = bases.choice()
			ref_base = base_pair(SHINE_DALGARNO[i])
			while rand_base == ref_base:
				rand_base = bases.choice()
			first_filler += rand_base

		first_filler = first_filler[::-1]

		riboswitch_sequence = first_filler + comp_base + revcomp_effector + revcomp_filler + loop + filler + dna_effector

		return riboswitch_sequence


	upstream_list = []

	#-35 sequence
	upstream_list.append(NEGATIVE_35_SEQUENCE)

	#find how long the filler needs to be, this will depend on the length of the -35 sequence
	index_end = -35 + len(NEGATIVE_35_SEQUENCE) #finds where -35 seq ends
	filler_length = -10 - index_end #finds length between end of -35 sequence and -10 sequence. this will be the filler length

	upstream_list.append(generate_random_seq(filler_length))

	upstream_list.append(PRIBNOW_BOX)

	#find how long the filler needs to be, this will depend on the length of the -10 sequence
	index_end = -10 + len(PRIBNOW_BOX) #finds where -10 seq ends
	filler_length = 0 - index_end #finds length between end of -10 sequence and Start Transcription Site

	upstream_list.append(generate_random_seq(filler_length))

	seq_list = parse_fasta(filename)
	effector = seq_list[0].sequence

	riboswitch_sequence = ''
	if 'a' not in effector:
		riboswitch_comps = generate_alternative_riboswitch(effector)
		riboswitch_comps = generate_riboswitch(effector)
		riboswitch_sequence = ''
		for components in riboswitch_comps:
			riboswitch_sequence += components

		LOOP_OFFSET = len(riboswitch_comps[0]) + len(riboswitch_comps[1]) - 1
		LOOP_LENGTH = len(riboswitch_comps[2])
			
		LENGTH_OF_X = len(riboswitch_sequence)


	else:
		riboswitch_comps = generate_riboswitch(effector)
		riboswitch_sequence = ''
		for components in riboswitch_comps:
			riboswitch_sequence += components

		LOOP_OFFSET = len(riboswitch_comps[0]) + len(riboswitch_comps[1]) - 1
		LOOP_LENGTH = len(riboswitch_comps[2])
			
		LENGTH_OF_X = len(riboswitch_sequence)


	upstream_list.append(riboswitch_sequence)

	upstream_list.append(SHINE_DALGARNO)

	index_end = -8 + len(SHINE_DALGARNO)
	filler_length = index_end - 0

	upstream_list.append(generate_random_seq(filler_length))
	upstream_sequence = ''
	upstream_sequence = ''.join(upstream_list)
	return upstream_sequence

def generate_random_seq(length, cg_rich=False):
	"""Generates a random DNA sequence with the given length. If cg_rich == True, then the random sequence
	will strongly favor in CG composition.

	length - length of sequence (int)
	cg_rich - determines if sequence will be CG rich or not (boolean)

	returns DNA sequence (string)
	"""
	sequence = ''
	if cg_rich:
		CG_num=length*3/4 # the number of C's and G's will always be at least 75 percent of the length of the entire sequence
		rand_CG_num=random.randint(CG_num,length) # the number of C's and G's will be between 75% and 100% of the entire sequence length
		CGpool=['C','G']
		for i in range(rand_CG_num):
			CGseq=random.choice(CGpool)
		
		AT_number=length-CG_num # the A's and T's will fill in the rest of the sequence
		ATpool=['A','T']
		for j in range(AT_number):
			ATseq=random.choice(ATpool)

		comb=CGseq+ATseq # combine CGseq and ATseq
		seq_list = []
		for c in comb:
			seq_list.append(c)
		sequence=''.join(random.sample(seq_list,len(seq_list))) # shuffle the order of nucleotide


	else :
		nuc_pool=['A','C','T','G']
		for i in range(length):
			sequence += random.choice(nuc_pool) # generate a sequence that does not have any restrction on the content ratio

	return sequence

def generate_random_loop(length):
	"""Generates a random DNA sequence that for the loop portion of a stemloop (hairpin loop).

	length - length of sequence (int)

	returns DNA sequence (string)
	"""

	"""Generates a random DNA sequence that for the loop portion of a stemloop (hairpin loop).

	length - length of sequence (int)

	returns DNA sequence (string)
	"""
	rand = random.randint(6,10) # randomly choose the length of the loop to be a number between 6 and 10
	# one half of the loop is generated first and then the other half is determined 
	# to make sure that when the segment forms a loop the sequence will not basepair itself
	half = rand/2
	first_half_seq = generate_random_seq(half)
	A_pool = ['C','G','T']
	C_pool = ['A','G','T']
	G_pool = ['A','C','T']
	T_pool = ['A','C','G']
	other_half=''
	for nuc in first_half_seq:
		if nuc == 'A':
			other_half += random.choice(A_pool)
		elif nuc == 'C':
			other_half += random.choice(C_pool)
		elif nuc == 'G':
			other_half += random.choice(G_pool)
		elif nuc == 'T':
			other_half += random.choice(T_pool)
	# because the entire sequence is to form a loop, the "other_half" needs to be reversed in its order
	second_half_seq = other_half[::-1] 

	# determine if rand is an odd number;
	# if it is an even number, concatenate "first_half_seq" and "second_half_seq" and output the resulting sequence;
	# if it is an odd number, add a randomly chosen nucleotide in between the "first_half_seq" and the "second_half_seq"
	# to generate an output sequence of the length, called "rand", which has been decided in the first line of the function
	if rand%2 == 0: 
		loop=first_half_seq+second_half_seq
	else :
		nuc_pool=['A','C','G','T']
		rand_nuc=random.choice(nuc_pool)
		loop=first_half_seq+rand_nuc+second_half_seq

	return loop



def generate_terminator():
	"""Generates the intrinsic termiantor sequence

	returns terminator sequence (string)
	"""

	# first produce a CG-rich segment (of a random length) and name it "first_half"
	# produce its reverse complement and name it "second_half";
	# generate a short sequence composed of randomly picked nucleotides (and of random length) to connect these two segments to form a hairpin loop and name it "middle";
	# generate a poly T tail of a random length and name it "poly_Ttail"; this will be complemented to produce poly U tail;
	# produce a random sequence to elongate one end of the terminator for the purpose of riboswitch (named "chunk");
	# concatenate the segments in the order of "chunk_seq"+first_half"+"middle"+"second_half"+"poly_Ttail"

	rand_first=random.randint(6,12) # randomly choose the length of the CG-rich portion of terminator, the sequence called "first_half"
	first_half = generate_random_seq(rand_first,True) # produce a CG-rich sequence 

	second_half=reverse_complement(first_half) # produce reverse complementary sequence of first_half 

	middle = generate_random_loop(random.randint(4,8)) # produce the loop portion 

	rand_Ttail=random.randint(8,12) # randomly choose the length of poly T tail
	poly_Ttail='T'*rand_Ttail # generate the poly T tail sequence

	rand_chunk=random.randint(15,20) # randomly choose the length of "chunk sequence"
	chunk_seq=generate_random_seq(rand_chunk) # generate the chunk sequence

	FILLER_SIZE = len(chunk_seq)

	terminator_sequence=chunk_seq+first_half+middle+second_half+poly_Ttail # concatenate the 5 segments to produce the finalized terminator sequence

	return terminator_sequence


def parse_params(filename):
	"""Given the file name of a parameters file, parse the file for parameter N.

	filename - file name of parameters file

	returns parameter N
	"""

	"""
	TO DO:
	-use parse fasta
	"""
	infile = None
	try:
		infile = open(filename, 'r')
	except IOError:
		print ("Error: Unable to open Parameters file, " + filename)
		sys.exit(1)

	#pulls the lines
	count = 0
	for line in infile:
		count += 1
		if count == 1:
			try:
				value = int(line)
			except ValueError:
				print("Error: Parameters not castable to int, " + line)
				sys.exit(1)
		elif (count >= 1 and line != None):
			raise IOError("Error: Parameter file not properly formed, " + line)
			sys.exit(1)
	return value



"""
MAIN FUNCTION
"""
"""
output = []

output.append(generate_upstream(EFFECTOR_FILENAME))

aa_seq = generate_aa_seq(PROTEIN_FILENAME)
coding_seq = convert_aa_to_dna(aa_seq, parse_codonfreq(CODON_FREQ_FILENAME))

output.append(coding_seq)

output.append(generate_terminator())
"""


#output = []

upstream_seq = generate_upstream(EFFECTOR_FILENAME)
upstreamlength = len(upstream_seq)
LENGTHS.append(upstreamlength)
#output.append(upstream_seq)



aa_seq = generate_aa_seq(parse_fasta(PROTEIN_FILENAME))
aaseqlength = len(aa_seq)
LENGTHS.append(aaseqlength)
coding_seq = convert_aa_to_dna(aa_seq, parse_codonfreq(CODON_FREQ_FILENAME))

#output.append(coding_seq)

term_seq = generate_terminator()
termseqlength = len(term_seq)
LENGTHS.append(termseqlength)
#output.append(term_seq)


N_variant = parse_params(PARAMS_FILENAME)

seq_list = parse_fasta(PROTEIN_FILENAME)

SITES = parse_sites(RESTRIC_ENZ_FILENAME)

sequence = generate_aa_seq(seq_list)
print(sequence)
sequence = convert_aa_to_dna(sequence, parse_codonfreq(CODON_FREQ_FILENAME))
print(sequence)

