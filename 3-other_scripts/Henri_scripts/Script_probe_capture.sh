

#  Script.sh
#  
#
#  Created by Henri Siljanen
#  

# Generation of hmmer profiles:
# The genes were collected with text based from NCBI nucletide database. The cultured organisms were also collected with blastN searches against nt-database without uncultured organims. The genes of cultured organims were collected, genes were aligned in the AliView program with MUSCLE alignment. These reference gene databases were used for In House formed hmmer-profile production.


# In Taito-cluster:
ssh taito-shell.csc.fi -l siljanen
cd $WRKDIR/

# Download the biokit which have the required tool:
# hmmer
# cd-hit


module load biokit

hmmbuild napA.hmm napA_refdb_NEW.fasta
hmmbuild narG.hmm narG_refdb_NEW.fasta
hmmbuild nirK.hmm nirK_refdb_NEW.fst
hmmbuild nosZ.hmm nosZ_refdb_NEW.fst
hmmbuild nxrB.hmm nxrB_refdb_NEW.fasta
hmmbuild TamoA.hmm TamoA_refdb_NEW.fst
hmmbuild nifH.hmm nifH_refdb_NEW.fst
hmmbuild nirS.hmm nirS_refdb_NEW.fst
hmmbuild nrfA.hmm nrfA_refdb_NEW.fa
hmmbuild pmoA.hmm pmoA_refdb_NEW.fasta
hmmbuild amoA.hmm amoA_bacterial_refDB_NEW.fst
hmmbuild hzoA.hmm hzoA_refdb_alignment.txt
hmmbuild mcrA.hmm mcrA_fungene_seed_aligned.fasta


# for norB gene, hmmer profile from fungene was used:
#norB.hmm
#mmoX.hmm
 
 
mv *.hmm $WRKDIR/DONOTREMOVE/compiled_seed



#################
# Search of text based from NCBI-nucleotide database, March 2017.
# files are downloaded to the desktop and compiled to files, text searched/hand-collected files:
napA_ncbi_new.fasta
pmoA_ncbi_new.fasta
narG_ncbi_new.fasta
TamoA_ncbi_new.fasta
nirk_ncbi_new.fasta
nirS_ncbi_new.fasta
nxrB_ncbi_new.fasta
amoA_ncbi_new.fasta
nrfA_ncbi_new.fasta
nifH_ncbi_new.fasta
hzoA_ncbi_new.fasta

# move databases to Taito
#folder- $WRKDIR/DONOTREMOVE/compiled_seed
# and search with hmmer from them.

#on Desktop
scp script_for-hmmer_NEW.sh siljanen@taito.csc.fi:/wrk/siljanen

ssh taito-shell.csc.fi -l siljanen
cd $WRKDIR/

#in Taito
sbatch script_for-hmmer_NEW.sh
#13/03/2017 klo 13:09
#Submitted batch job 12800818


# contect of: script_for-hmmer_NEW.sh

#!/bin/bash -l
#SBATCH -J hmmer_against_ncbi_hand-collected_job
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 72:00:00
#SBATCH -n 4
#SBATCH --mail-type=END
#SBATCH --mail-user=henri.siljanen@uef.fi
#SBATCH --mem-per-cpu=8000
#


cd $WRKDIR/DONOTREMOVE/compiled_seed

nhmmer --noali --cpu 4 napA.hmm napA_ncbi_new.fasta > napA_hmmer_list3.txt
nhmmer --noali --cpu 4 pmoA.hmm pmoA_ncbi_new.fasta > pmoA_hmmer_list3.txt
nhmmer --noali --cpu 4 narG.hmm narG_ncbi_new.fasta > narG_hmmer_list3.txt
nhmmer --noali --cpu 4 TamoA.hmm TamoA_ncbi_new.fasta > TamoA_hmmer_list3_2.txt
nhmmer --noali --cpu 4 nosZ.hmm nosZ_ncbi_new.fasta > nosZ_hmmer_list3.txt
nhmmer --noali --cpu 4 nirK.hmm nirk_ncbi_new.fasta > nirK_hmmer_list3.txt
nhmmer --noali --cpu 4 nirS.hmm nirS_ncbi_new.fasta > nirS_hmmer_list3.txt
nhmmer --noali --cpu 4 nxrB.hmm nxrB_ncbi_new.fasta > nxrB_hmmer_list3.txt
nhmmer --noali --cpu 4 amoA.hmm amoA_ncbi_new.fasta > amoA_hmmer_list3.txt
nhmmer --noali --cpu 4 nrfA.hmm nrfA_ncbi_new.fasta > nrfA_hmmer_list3.txt
nhmmer --noali --cpu 4 nifH.hmm nifH_ncbi_new.fasta > nifH_hmmer_list3.txt

nhmmer --noali --cpu 4 hzoA.hmm hzoA_ncbi_new.fasta > TamoA_hmmer_list3.txt
nhmmer --noali --cpu 4 mcrA.hmm mcrA_ncbi_new.fasta > mcrA_hmmer_list3.txt


#################################################################
# Fungene databases download. Genes downloaded from Fungene 8.5.
# in Downloads
#
# example compilement of napA with cat-command

cat fungene_8.5_napA_10000_unaligned_nucleotide_seqs-2.fa fungene_8.5_napA_10000_unaligned_nucleotide_seqs-3.fa fungene_8.5_napA_10000_unaligned_nucleotide_seqs.fa fungene_8.5_napA_511_unaligned_nucleotide_seqs.fa fungene_8.5_napA_9021_unaligned_nucleotide_seqs.fa > napA_all_fungene.txt

# Transfer of file to Taito.
scp napA_all_fungene.txt siljanen@taito.csc.fi:/wrk/siljanen/DONOTREMOVE/compiled_seed


# July 2017:

#now new seed files and hmm-profiles are in  $WRKDIR/DONOTREMOVE/compiled_seed

ssh taito-shell.csc.fi -l siljanen


cd $WRKDIR/DONOTREMOVE/compiled_seed

# get the nt and nr database in Taito

wget  ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz




###############################################

#in Desktop of MacBook pro/ transfer .sh file to Taito-cluster
scp script_for_hmmer_against_nt.sh siljanen@taito.csc.fi:/wrk/siljanen

# contect of: script_for_hmmer_against_nt.sh

#!/bin/bash -l
#SBATCH -J hmmer_against_NT_job
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 72:00:00
#SBATCH -n 4
#SBATCH --mail-type=END
#SBATCH --mail-user=henri.siljanen@uef.fi
#SBATCH --mem-per-cpu=8000
#


cd $WRKDIR/DONOTREMOVE/compiled_seed

gunzip nt.gz

module load biokit

nhmmer --noali --cpu 4 hzoA.hmm nt > hzoA_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 amoA.hmm nt > amoA_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 mcrA.hmm nt > mcrA_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 napA.hmm nt > napA_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 nifH.hmm nt > nifHA_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 nirS.hmm nt > nirS_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 nrfA.hmm nt > nrfA_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 pmoA.hmm nt > pmoA_hmmer_results_from_nt.txt

nhmmer --noali --cpu 4 mmoX-trimmed.hmm nt > mmoX_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 narG.hmm nt > narG_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 nirK.hmm nt > nirK_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 nosZ.hmm nt > nosZ_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 nxrB.hmm nt > nxrB_hmmer_results_from_nt.txt
nhmmer --noali --cpu 4 TamoA.hmm nt > TamoA_hmmer_results_from_nt.txt


##############################################

# in CSC Taito server
# in Taito  --- submit blastn run from nt database:


cd $WRKDIR
script_for_hmmer_against_nt.sh


sbatch script_for_hmmer_against_nt.sh
Submitted batch job 15339474 27/07-17 at 00:37am



###############################################

# download env_nt database of all metagenomic projects  of WGS (previous SRA) database search:
# in Taito: folder $WRKDIR/DONOTREMOVE/compiled_seed

wget  ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/env_nt.gz

gunzip env_nt.gz


###############


#In MacBook pro:
scp script_for_hmmer_against_env_nt.sh siljanen@taito.csc.fi:/wrk/siljanen

#Content of: script_for_hmmer_against_env_nt.sh

#!/bin/bash -l
#SBATCH -J hmmer_against_env_nt_job
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 72:00:00
#SBATCH -n 4
#SBATCH --mail-type=END
#SBATCH --mail-user=henri.siljanen@uef.fi
#SBATCH --mem-per-cpu=8000
#

cd $WRKDIR/DONOTREMOVE/compiled_seed

module load biokit

nhmmer --noali --cpu 4 hzoA.hmm env_nt > hzoA_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 amoA.hmm env_nt > hzoA_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 mcrA.hmm env_nt > mcrA_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 napA.hmm env_nt > napA_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 nifH.hmm env_nt > nifH_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 nirS.hmm env_nt > nirS_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 nrfA.hmm env_nt > nrfA_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 pmoA.hmm env_nt > pmoA_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 norB.hmm env_nt > norB_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 mmoX.hmm env_nt > mmoX_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 narG.hmm env_nt > narG_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 nirK.hmm env_nt > nirK_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 nosZ.hmm env_nt > nosZ_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 nxrB.hmm env_nt > nxrB_hmmer_results_from_env_nt.txt
nhmmer --noali --cpu 4 TamoA.hmm env_nt > TamoA_hmmer_results_from_env_nt.txt


# in CSC Taito server
# in Taito  --- submit blastn run from env_nt database:

cd $WRKDIR

sbatch script_for_hmmer_against_env_nt.sh
#Submitted batch job 15339472 27/07-17 at 00:37am


### COMPILING all files together:
# I copied seed data after fungene compilement and my own hmmer runs and also blast searches from wgs databases
# into just one folder $WRKDIR/DONOTREMOVE/compiled_seed

# Henri's hmmer files from NCBI nucleotide, nt and env_nt databases, and fungene compiled wth cat command

cat amoA_hmmer_list3.txt amoA_hmmer_results_from_nt.txt amoA_hmmer_results_from_env_nt.txt amoA_all_fungene.txt > amoA_ncbi_nt_env_nt_fungene_all_compiled.txt

cat TamoA_hmmer_list3.txt TamoA_hmmer_results_from_nt.txt TamoA_hmmer_results_from_env_nt.txt TamoA_all_fungene.txt > TamoA_ncbi_nt_env_nt_fungene_all_compiled.txt

cat nxrB_hmmer_list3.txt nxrB_hmmer_results_from_nt.txt nxrB_hmmer_results_from_env_nt.txt nxrB_all_fungene.txt > nxrB_ncbi_nt_env_nt_fungene_all_compiled.txt

cat napA_hmmer_list3.txt napA_hmmer_results_from_nt.txt napA_hmmer_results_from_env_nt.txt napA_all_fungene.txt > napA_ncbi_nt_env_nt_fungene_all_compiled.txt

cat narG_hmmer_list3.txt narG_hmmer_results_from_nt.txt narG_hmmer_results_from_env_nt.txt narG_all_fungene.txt > narG_ncbi_nt_env_nt_fungene_all_compiled.txt

cat nirS_hmmer_list3.txt nirS_hmmer_results_from_nt.txt nirS_hmmer_results_from_env_nt.txt nirS_all_fungene.txt > nirS_ncbi_nt_env_nt_fungene_all_compiled.txt
cat nirK_hmmer_list3.txt nirK_hmmer_results_from_nt.txt nirK_hmmer_results_from_env_nt.txt nirK_all_fungene.txt > nirK_ncbi_nt_env_nt_fungene_all_compiled.txt

cat norB_hmmer_list3.txt norB_hmmer_results_from_nt.txt norB_hmmer_results_from_env_nt.txt norB_all_fungene.txt > norB_ncbi_nt_env_nt_fungene_all_compiled.txt
cat nosZ_hmmer_list3.txt nosZ_hmmer_results_from_nt.txt nosZ_hmmer_results_from_env_nt.txt nosZ_all_fungene.txt nosZ_nuc_derep_unmasked.fasta > nosZ_ncbi_nt_env_nt_fungene_S-Hallin_all_compiled.txt
cat nifH_hmmer_list3.txt nifH_hmmer_results_from_nt.txt nifH_hmmer_results_from_env_nt.txt nifH_all_fungene.txt > nifH_ncbi_nt_env_nt_fungene_all_compiled.txt
cat pmoA_hmmer_list3.txt pmoA_hmmer_results_from_nt.txt pmoA_hmmer_results_from_env_nt.txt pmoA_all_fungene.txt > pmoA_ncbi_nt_env_nt_fungene_all_compiled.txt
cat mmoX_hmmer_list3.txt mmoX_hmmer_results_from_nt.txt mmoX_hmmer_results_from_env_nt.txt mmoX_all_fungene.txt > mmoX_ncbi_nt_env_nt_fungene_all_compiled.txt
cat mcrA_hmmer_list3.txt mcrA_hmmer_results_from_nt.txt mcrA_hmmer_results_from_env_nt.txt mcrA_all_fungene.txt > mcrA_ncbi_nt_env_nt_fungene_all_compiled.txt

cat hzoA_hmmer_list3.txt hzoA_hmmer_results_from_nt.txt hozA_hmmer_results_from_env_nt.txt hzoA_all_fungene.txt > hzoA_ncbi_nt_env_nt_fungene_all_compiled.txt

#compile all the dataset together:

cat *_all_compiled.txt > all_genes_compiled_ncbi_nt_envnt_fungene.fasta

#cluster sequences to unique ones:

cd-hit-est -i all_genes_compiled_ncbi_nt_envnt_fungene.fasta -o  all_genes_compiled_ncbi_nt_envnt_fungene_unique.fasta -c 1.00

# SeqCap probe design:

# all_genes_compiled_ncbi_nt_envnt_fungene_unique.fasta , were submitted to SeqCap website, with 80% clustering and 6 probes per cluster:
 
 
#####################################################

# Angus:

# Check up of amoA vs. pmoA reads. Also whether the narG or napA reads have formate dehydrogenase. This then produced corrected dataframe for the plots. I think you had the script for doing the comparison with blast and then rejecting those reads which were not amoA, or pmoA, or narG, or napA.

#The TF/FP calculations:


###############################################################################
# Statitical test for Origical and Targeted reads (Mock community comparison) , plotting the Figure (Script_Fig_S2.sh / Script_ProbeCapture_Fig_S3.sh)
# Statitical test for Shotgun and Targeted reads (Agricultual and wetland soil), plotting the Figure  (Script_Fig_2.sh / Script_Fig_3.sh).



