#for parsing HMMer 3.0 hmmscan table output
import string
import argparse
import os
import re
import subprocess
from Bio.Seq import Seq
from Bio import SeqFeature
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIWWW
from Bio import Entrez
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from urllib2 import HTTPError
from urllib2 import URLError

import random

Entrez.email = 'aurelien.saghai@slu.se'
Entrez.api_key = '2549482f72f61af0de0df3cacb06385e8e09'

def translate_six_frames(seq):
	"""Translate nucleotide sequences to all six frames"""
	
	rev = seq.reverse_complement()
	for i in range(3):
		yield seq[i:]
		yield rev[i:]

def convert_seq_to_record(seq,seqid):
	"""Converts sequence to nucleotide and amino acid sequences for HMMER"""
	
	nuc_frame = SeqRecord(seq)
	prot_frame = SeqRecord(seq.translate(table=1))
	nuc_frame.id = seqid
	prot_frame.id = seqid
	return (nuc_frame,prot_frame)

def split_len(seq, length):
	return [seq[i:i+length] for i in range(0, len(seq), length)]

def rand_tag():
	charList = string.letters+string.digits
	tag=''
	for i in range(0,5):
		tag += str(charList[random.randint(0,len(charList)-1)])
	return tag
		

#get user files and arguments
parser = argparse.ArgumentParser(
description = '''Script for aligning nucleotide sequences using a seed amino acid alignment and HMMER.''')
parser.add_argument('-f',help='FASTA File of nucleotide sequences to be aligned.')
parser.add_argument('-s',help='FASTA file with seed amino acid alignment. The HMM is generated from this file')
parser.add_argument('-o',help='FASTA file for output alignment.')
parser.add_argument('-e',help='e-value cutoff. If input sequence has no domains that match the HMM at this e-value threshold, it is placed in the error_seq.fasta file for inspection.')
parser.add_argument('-d',help='Dereplicate amino acid sequences at <d> similarity. Must be value between 0 and 1. If left blank (default) no dereplication is performed and all sequences are kept')
args = vars(parser.parse_args())

seed = os.getcwd()+'/'+args['s']
inseq_file = os.getcwd()+'/'+args['f']
evalue = args['e']
out_file = os.getcwd()+'/'+args['o']
derep = args['d']

#initiating lists for objects
prot = []
nuc = []
seq_accs = []
seq_error=[]
seq_desc={}

############################################################################
#generate HMM model for scanning and alignment
###########################################################################
seed_in = SeqIO.parse(seed, "fasta")
out_handle = open("hmm_seed.sto","w")
SeqIO.write(seed_in, out_handle, 'stockholm')
out_handle.close()

print('###########################################################################\nGenerating HMM model\n###########################################################################')
os.system('hmmbuild hmm.out hmm_seed.sto')
os.system('hmmpress hmm.out')
############################################################################
#screen input sequences for frameshift errors based on AA search, 6-frame
###########################################################################

print('###########################################################################\nScreening sequences with e-value of:'+evalue)

error_log = open('error_log.txt', 'w')
frame_log = open('framshift_log.txt', 'w')

nseq = 0
count = 0
for seq_record in SeqIO.parse(inseq_file, "fasta"):

	#trimming garbage off of assembly accession name
	check_sid = seq_record.id.split('_')
	#fix to ensure dots are removed properly
	if len(check_sid) == 1:
	
		sid = check_sid[0]
	
	else:
		if re.search('\.',check_sid[1]):
	
			if len(check_sid) > 1:
				sid = check_sid[0]+'_'+check_sid[1]
			else:
				sid = check_sid[0]

		else:
			#fix to snure dashes are kept.
			sid = seq_record.id
			#print sid
	
	seq_record.id = sid
	
	if seq_record.id in seq_accs:
		
		#add integer to assembly accession name to indicate multiple gene instances in a genome
		res = len([i for i in seq_accs if seq_record.id in i])+1
		new = seq_record.id+'_'+str(res)
		error_log.write(seq_record.id+' - duplicated id, changed to '+new+'\n')
		seq_record.id = new

	seq_accs.append(seq_record.id)
	seq_record.seq =  Seq(re.sub('[\*-?\.]','',str(seq_record.seq)))
	
	try:
		seq_record.seq.translate(table=1)
	except:
		error_log.write(seq_record.id+" skipped - Could not translate after gap removal"+'\n')
		continue

	seq_out = []
	#translate all six frames
	frame = list(translate_six_frames(seq_record.seq))
	#convert each frame to seq object

	fi=1
	for f in frame:	
		#if stop codon in middle next frame
			
			s_out = SeqRecord(f.translate(table=1))
			s_out.id = 'frame'+str(fi)
			seq_out.append(s_out)
			fi=fi+1
			#true frame is candidate with highest HMM score

	#write all frames to fasta for HMMER - run HMMscan to get probable domains...
	SeqIO.write(seq_out,"pytest.fasta","fasta")
	z = subprocess.Popen('hmmscan --domtblout dom.out -E '+evalue+' --domE '+evalue+' -Z 1000000 hmm.out pytest.fasta',shell=True, stdout=subprocess.PIPE).communicate()
	file = open('dom.out','r+')

	#get start and end coordinates as well as frames from hmmer output
	starts = []
	ends = []
	frames = []

	for line in file:
		if ('#') in line:
			continue
		else:
			#the coordinate scores are based on the 'ali pos' in the HMMscan output
			coord = re.sub('\s+','|',str(line)).split('|')

			#need to be sure to subtract 1 for positions in python
			start = (int(coord[17])-1)*3
			end = (int(coord[18])-1)*3
			fr = int(re.sub('frame','',coord[3]))-1
			frseq = str(list(frame)[int(re.sub('frame','',coord[3]))-1])
			starts.append(start)
			ends.append(end)
			frames.append(fr)
	#if no significant frames, send to error log
	if len(frames) == 0:
		error_log.write(seq_record.id+' - '+seq_record.description+' - No domains below E = '+evalue+'\n')
		seq_error.append(seq_record)
		continue
	#if many error frames, send to error log
	elif len(frames) > 10:
		error_log.write(seq_record.id+' - '+seq_record.description+' - Many frameshift errors'+'\n')
		continue

	#if only one reading frame, write nuc to one file and translation to another
	elif len(frames) == 1:
		
		out = convert_seq_to_record(frame[fr],seq_record.id)
		
		prot.append(out[1])
		nuc.append(out[0])

	#otherwise for seqs with several frames, concatenate sequences and signify framshift location
	#generate breakpoints
	else:
		
		#organize all frames and breakpints in the proper order
		fr_seg = []
		breaks = []
		#find first frame based on starting point in alignment 
		ind = starts.index(min(starts))
		breaks.append(0)
		breaks.append(ends[ind])
		fr_seg.append(frames[ind])
		
		#pop first fragment, iterate through the rest
		starts.pop(ind)
		ends.pop(ind)		
		frames.pop(ind)

		while len(frames) > 0:
			#find subsequent frame and add to list
			ind = starts.index(min(starts))
			indb = breaks.index(max(breaks))
			
			#if start of this fragment is before the end of the previous
			if breaks[indb] > starts[ind]:
				
				#make end of previous fragment the same coordinates as the start of this fragment
				breaks[indb] = starts[ind]
				breaks.append(ends[ind])
			else:
				breaks.append(ends[ind])
				
			fr_seg.append(frames[ind])
			starts.pop(ind)
			ends.pop(ind)
			frames.pop(ind)

		segs = []
		
		frame_log.write(seq_record.id+' - Frameshift detected in positions '+str(breaks)+'\n')
		for b in range(1,len(breaks)):
			#get each fragment based on frame and break points
			segs.append(str(frame[fr_seg[b-1]][breaks[b-1]+3:breaks[b]]))
		
		#join fragments and add 'NNN' to indicate join in alignment
		out = convert_seq_to_record(Seq("NNN".join(segs)),seq_record.id)
		
		nuc.append(out[0])
		prot.append(out[1])
		
	nseq += 1
	seq_desc[seq_record.id] = seq_record.description
#create dictionary to call sequences by key for nucleotide mapping
nuc_dict = SeqIO.to_dict(nuc)

###################################################################
#then align protein by hmmer and map to nucleotide in fasta format
prot_handle = open("protout.fasta", "w")
SeqIO.write(prot,prot_handle,"fasta")

prot_handle.close()
SeqIO.write(seq_error,"error_seq.fasta","fasta")
print('\nAligning '+str(nseq)+' amino acid sequences by HMM')

if derep != None:

	print('\nRunning CD-HIT at '+str(derep)+' amino acid sequences by HMM')
	subprocess.Popen('cd-hit -T 0 -c ' +str(derep)+ ' -i ./protout.fasta -o ./protout_derep.fasta',shell=True).communicate()
	subprocess.Popen('hmmalign -o hmmaligned.sto hmm.out ./protout_derep.fasta',shell=True, stdout=subprocess.PIPE).communicate()
	

else:

	subprocess.Popen('hmmalign -o hmmaligned.sto hmm.out protout.fasta',shell=True, stdout=subprocess.PIPE).communicate()


ali_in = AlignIO.read(open("hmmaligned.sto"), "stockholm")
out_handle = open("prot_aligned.fasta","w")
AlignIO.write(ali_in, out_handle, 'fasta')
out_handle.close()

nuc_aligned = []
gb_aligned = []
print('###########################################################################\nMapping nucleotides to aligned amino acid sequences')
###################################################################
#now assign if user has specified
	
for seq in ali_in:

	prot_seq =  seq.seq
	
	nuc_seq = nuc_dict[seq.id].seq
	pr = nuc_seq.translate(table=1)

	codons = split_len(str(nuc_seq),3)
	codon_aligned = []
	pos = 0
	
	for i in range(pos,len(prot_seq)):
		
		if prot_seq[i] == '-':
			codon_aligned.append('---')
		elif prot_seq[i] == '.':
			codon_aligned.append('...')
		
		else:
			
			aa = Seq(codons[pos]).translate(table=1,stop_symbol='X')
			#print aa
			aa_check = pr[pos]
			
			if prot_seq[i].islower():

				codon_aligned.append(codons[pos].lower())

			else:
				codon_aligned.append(codons[pos].upper())
			pos += 1

	nuc = Seq(''.join(codon_aligned))
	nuc.alphabet = generic_dna	
	seq_record = nuc_dict[seq.id]
	seq_record.seq = nuc
	seq_record.description = seq_desc[seq.id]

	nuc_aligned.append(seq_record)
	
nuca_handle = open(out_file, "w")
SeqIO.write(nuc_aligned,nuca_handle,"fasta")
nuca_handle.close()
error_log.close()
frame_log.close()

subprocess.Popen('rm hmm_seed.sto hmm.out* hmmaligned.sto protout.fasta pytest.fasta dom.out protout_derep.fasta.clstr',shell=True, stdout=subprocess.PIPE).communicate()

