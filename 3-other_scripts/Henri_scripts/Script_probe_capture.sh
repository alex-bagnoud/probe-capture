

#  Script.sh
#  
#
#  Created by Henri Siljanen
#  

# Generation of hmmer profiles:
# 1°. The genes were collected with text based from NCBI nucletide database. 
#         1.1°. Search prases were used to search the sequences from different phylogenetic groups of full size of each genes. 
#          1.2°. The cultured organisms were also collected with blastN searches against nt-database without uncultured organims. 
#                   This database was used for production of hmmer-profile. 
# 2°. The genes of cultured organims were aligned in the AliView program with MUSCLE alignment. 
#          These reference gene databases were used for In House formed hmmer-profile production.
# 3°. Search of each genes. The text based search from NCBI produces also full genomic fragments, in step 1°. 
#            The produced House formed hmmer profile were used to search each genes from these text based searches in . 
# 4°. The fungene database for each genes were searched. 
# 5°. The full nt and env-nt database was downladed on the server, from NCBI. These database searched with House made hmmer-profiles to build the database 
#      5.1°.from nt 
#      5.2°. and env-nt.  
# 6°. Compiling files: These files from NCBI text search, fungene, private databse and hmmer searches from nt and env-nt databse were compiled together 
#         and submitted to MetCap pipeline to produce probes with Roche SeqCap method.  
# 7°. Submitting Target gene database to MetCap probe design:
# 8°. Submitting probes to SeqCap probe order:
#

#################
# 1°. Search of text based from NCBI-nucleotide database, March 2017.
# Text based seacrh terms for each gene, NCBI website https://www.ncbi.nlm.nih.gov/nuccore : 

#mmoX 
"soluble methane monooxygenase" 
mmoX[gene name] NOT "soluble methane monooxygenase"
#pmoA
"particulate methane monooxygenase subunit A"
pmoA[gene name] NOT "particulate methane monooxygenase subunit A"
#amoA 
"ammonia monooxygenase subunit A"
amoA[gene name] NOT "ammonia monooxygenase subunit A"
# narG
“nitrate reductase alpha chain”
“nitrate reductase alpha subunit”
narG[gene name]
“nitrate reductase alpha chain” NOT “nitrate reductase alpha subunit” NOT narG[gene name]
“nitrate reductase catalytic subunit” NOT “nitrate reductase alpha chain” NOT narG[gene name]
“nitrate reductase alpha subunit” NOT narG[gene name] 
#napA
"periplasmic nitrate reductase"
napA[gene name]
#nifH
“nitrogenase reductase”
nifH[gene name]
“nitrogenase reductase” NOT nifH[gene name]
#nirK
“copper-dependent nitrite reductase”
nirK[gene name]
“copper-dependent nitrite reductase” NOT nirK[gene name]
#nirS
“cytochrome cd1 nitrite reductase”
nirS[gene name]
“cytochrome cd1 nitrite reductase” NOT nirS[gene name]
#nosZ
nosZ[gene name]
“nitrous oxide reductase” NOT nosZ[gene name]
putative nitrous oxide reductase NOT nosZ[gene name] 
#nrfA
nrfA[gene name]
“cytochrome c nitrite reductase” NOT nrfA[gene name]
#nxrB
nxrB[gene name]
nxrB[gene name] NOT "nitrite oxidoreductase beta subunit”
#hzoA
hydrazine-oxidizing enzyme NOT hzo[gene name]
hydrazine oxidoreductase NOT hzo[gene name]
hzo[gene name]
#mcrA
mcrA[gene name]
"methyl coenzyme M reductase"
"methyl coenzyme M reductase" NOT mcrA[gene name]

# Sequences files were downloaded to the desktop and compiled to files, text searched/hand-collected files:
# These searches produced several sequence files. As example mcrA produced three different files, those were merged together with cat-command: 
cat mcrA_ncbi_search_part1.txt  mcrA_ncbi_search_part2.txt  mcrA_ncbi_search_part2.txt > mcrA_ncbi_new.fasta

nifH_ncbi_new.fasta
TamoA_ncbi_new.fasta
amoA_ncbi_new.fasta
nxrB_ncbi_new.fasta
hzoA_ncbi_new.fasta
nrfA_ncbi_new.fasta
napA_ncbi_new.fasta
narG_ncbi_new.fasta
nirK_ncbi_new.fasta
nirS_ncbi_new.fasta
nosZ_ncbi_new.fasta
pmoA_ncbi_new.fasta
mmoX_ncbi_new.fasta
mcrA_ncbi_new.fasta

# 1.1.° During these searches, the results of each taxon was collected for separated search file. 

# 1.2.°  These sequences were used to search similar cultured organisms in BlastN search (exclution of uncultured sequences). 
# https://blast.ncbi.nlm.nih.gov/ 

# 2°.This search produced the simple databases for House made hmmer profiles. The sequences were aligned with AliView programs with MUSCLE algoritm. 

nifH_refdb_NEW.fst
TamoA_refdb_NEW.fst
amoA_bacterial_refDB_NEW.fst 
nxrB_refdb_NEW.fasta
hzoA_refdb_alignment.txt
nrfA_refdb_NEW.fa
napA_refdb_NEW.fasta
narG_refdb_NEW.fasta
nirK_refdb_NEW.fst
nirS_refdb_NEW.fst
nosZ_refdb_NEW.fst
pmoA_refdb_NEW.fasta

# In Taito-cluster:
ssh taito-shell.csc.fi -l siljanen
cd $WRKDIR/

# Download the biokit which have the required tool:
# hmmer
# cd-hit

module load biokit

hmmbuild nifH.hmm nifH_refdb_NEW.fst
hmmbuild TamoA.hmm TamoA_refdb_NEW.fst
hmmbuild amoA.hmm amoA_bacterial_refDB_NEW.fst
hmmbuild nxrB.hmm nxrB_refdb_NEW.fasta
hmmbuild hzoA.hmm hzoA_refdb_alignment.txt
hmmbuild nrfA.hmm nrfA_refdb_NEW.fa
hmmbuild napA.hmm napA_refdb_NEW.fasta
hmmbuild narG.hmm narG_refdb_NEW.fasta
hmmbuild nirK.hmm nirK_refdb_NEW.fst
hmmbuild nirS.hmm nirS_refdb_NEW.fst
hmmbuild nosZ.hmm nosZ_refdb_NEW.fst
hmmbuild pmoA.hmm pmoA_refdb_NEW.fasta

# for norB and mmoX gene, hmmer profile from fungene was used:
# norB.hmm
# mmoX.hmm 
# the fungene seed database for mcrA was used 
hmmbuild mcrA.hmm mcrA_fungene_seed_aligned.fasta
 
mv *.hmm $WRKDIR/DONOTREMOVE/compiled_seed


# 3°. Search of each genes. The text based search from NCBI produces also full genomic fragments, in step 1°. 
#            The produced House formed hmmer profile were used to search each genes from these text based searches in . 

# Move NCBI text based databases (nifH_ncbi_new.fasta , etc.) to Taito-supercomputer server
# folder- $WRKDIR/DONOTREMOVE/compiled_seed
# and search with hmmer the each genes from text search database.

#on Desktop
scp script_for_hmmer_NEW.sh siljanen@taito.csc.fi:/wrk/siljanen

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

# The nhmmer search of each gene from NCBI text based searches. Default values of nhmmer use (incE is 0.01). 
nhmmer --noali --cpu 4 nifH.hmm nifH_ncbi_new.fasta > nifH_hmmer_list3.txt
nhmmer --noali --cpu 4 amoA.hmm amoA_ncbi_new.fasta > amoA_hmmer_list3.txt
nhmmer --noali --cpu 4 TamoA.hmm TamoA_ncbi_new.fasta > TamoA_hmmer_list3_2.txt
nhmmer --noali --cpu 4 nxrB.hmm nxrB_ncbi_new.fasta > nxrB_hmmer_list3.txt
nhmmer --noali --cpu 4 hzoA.hmm hzoA_ncbi_new.fasta > hzoA_hmmer_list3.txt
nhmmer --noali --cpu 4 nrfA.hmm nrfA_ncbi_new.fasta > nrfA_hmmer_list3.txt
nhmmer --noali --cpu 4 napA.hmm napA_ncbi_new.fasta > napA_hmmer_list3.txt
nhmmer --noali --cpu 4 narG.hmm narG_ncbi_new.fasta > narG_hmmer_list3.txt
nhmmer --noali --cpu 4 nirK.hmm nirk_ncbi_new.fasta > nirK_hmmer_list3.txt
nhmmer --noali --cpu 4 nirS.hmm nirS_ncbi_new.fasta > nirS_hmmer_list3.txt
nhmmer --noali --cpu 4 nosZ.hmm nosZ_ncbi_new.fasta > nosZ_hmmer_list3.txt

nhmmer --noali --cpu 4 pmoA.hmm pmoA_ncbi_new.fasta > pmoA_hmmer_list3.txt
nhmmer --noali --cpu 4 mmoX.hmm mmoX_ncbi_new.fasta > mmoX_hmmer_list3.txt
nhmmer --noali --cpu 4 mcrA.hmm mcrA_ncbi_new.fasta > mcrA_hmmer_list3.txt


#################################################################
# 4°. The fungene database for each genes were searched.
# Fungene databases download. Genes downloaded from Fungene 8.5.
# in Downloads
#
# example compilement of napA with cat-command

cat fungene_8.5_napA_10000_unaligned_nucleotide_seqs-2.fa fungene_8.5_napA_10000_unaligned_nucleotide_seqs-3.fa fungene_8.5_napA_10000_unaligned_nucleotide_seqs.fa fungene_8.5_napA_511_unaligned_nucleotide_seqs.fa fungene_8.5_napA_9021_unaligned_nucleotide_seqs.fa > napA_all_fungene.txt

 # similarly files for 
 # nifH_all_fungene.txt
 # amoA_all_fungene.txt
 # TamoA_all_fungene.txt
 # nxrB_all_fungene.txt
 # nrfA_all_fungene.txt
 # nirK_all_fungene.txt
 # nirS_all_fungene.txt
 # norB_all_fungene.txt
 # nosZ_all_fungene.txt
 # mmoX_all_fungene.txt
 # pmoA_all_fungene.txt
 # mcrA_all_fungene.txt
 

# Transfer of file to Taito.
scp napA_all_fungene.txt siljanen@taito.csc.fi:/wrk/siljanen/DONOTREMOVE/compiled_seed
# same for other genes too. 

# 5°. The full nt and env-nt database was downladed on the server, from NCBI. These database searched with House made hmmer-profiles to build the database 
#      5.1°.from nt 
#      5.2°. and env-nt.  
# Search done on July 2017:

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
# 6°. Compiling files: These files from NCBI text search, fungene, private databse and hmmer searches from nt and env-nt databse were compiled together 
#         and submitted to MetCap pipeline to produce probes with Roche SeqCap method. 

# nifH database: ARB database from Zehr’s lab: http://wwwzehr.pmc.ucsc.edu/nifH_Database_Public/. 
# nifH_2014April04-2_export.fasta

# archaeal amoA: private database were compiled into the database 
#  Arch_amoA_clust96_20141111.txt (provider Ricardo Alves, emailed)

# pmoA database from Susanna Liebner https://dataservices.gfz-potsdam.de/panmetaworks/showshort.php?id=escidoc:1423157
# pmoa7809.fasta  
# I copied seed data after fungene compilement and my own hmmer runs and also blast searches from wgs databases
# into just one folder $WRKDIR/DONOTREMOVE/compiled_seed

# nosZ database: Private database from Sara Hallin and Chris Jones (get over email) nosZ_nuc_derep_unmasked.fasta

# Henri's hmmer files from NCBI nucleotide, nt and env_nt databases, and fungene compiled wth cat command

cat nifH_hmmer_list3.txt nifH_hmmer_results_from_nt.txt nifH_hmmer_results_from_env_nt.txt nifH_all_fungene.txt nifH_2014April04-2_export.fasta > nifH_ncbi_nt_env_nt_fungene_Zehr_all_compiled.txt
cat amoA_hmmer_list3.txt amoA_hmmer_results_from_nt.txt amoA_hmmer_results_from_env_nt.txt amoA_all_fungene.txt > amoA_ncbi_nt_env_nt_fungene_all_compiled.txt
cat TamoA_hmmer_list3.txt TamoA_hmmer_results_from_nt.txt TamoA_hmmer_results_from_env_nt.txt TamoA_all_fungene.txt Arch_amoA_clust96_20141111.txt > TamoA_ncbi_nt_env_nt_fungene_Alves_all_compiled.txt
cat nxrB_hmmer_list3.txt nxrB_hmmer_results_from_nt.txt nxrB_hmmer_results_from_env_nt.txt nxrB_all_fungene.txt > nxrB_ncbi_nt_env_nt_fungene_all_compiled.txt
cat hzoA_hmmer_list3.txt hzoA_hmmer_results_from_nt.txt hozA_hmmer_results_from_env_nt.txt hzoA_all_fungene.txt > hzoA_ncbi_nt_env_nt_fungene_all_compiled.txt
cat nrfA_hmmer_list3.txt nrfA_hmmer_results_from_nt.txt nrfA_hmmer_results_from_env_nt.txt nrfA_all_fungene.txt > nrfA_ncbi_nt_env_nt_fungene_all_compiled.txt
cat napA_hmmer_list3.txt napA_hmmer_results_from_nt.txt napA_hmmer_results_from_env_nt.txt napA_all_fungene.txt > napA_ncbi_nt_env_nt_fungene_all_compiled.txt
cat narG_hmmer_list3.txt narG_hmmer_results_from_nt.txt narG_hmmer_results_from_env_nt.txt narG_all_fungene.txt > narG_ncbi_nt_env_nt_fungene_all_compiled.txt
cat nirS_hmmer_list3.txt nirS_hmmer_results_from_nt.txt nirS_hmmer_results_from_env_nt.txt nirS_all_fungene.txt > nirS_ncbi_nt_env_nt_fungene_all_compiled.txt
cat nirK_hmmer_list3.txt nirK_hmmer_results_from_nt.txt nirK_hmmer_results_from_env_nt.txt nirK_all_fungene.txt > nirK_ncbi_nt_env_nt_fungene_all_compiled.txt
cat norB_hmmer_list3.txt norB_hmmer_results_from_nt.txt norB_hmmer_results_from_env_nt.txt norB_all_fungene.txt > norB_ncbi_nt_env_nt_fungene_all_compiled.txt
cat nosZ_hmmer_list3.txt nosZ_hmmer_results_from_nt.txt nosZ_hmmer_results_from_env_nt.txt nosZ_all_fungene.txt nosZ_nuc_derep_unmasked.fasta > nosZ_ncbi_nt_env_nt_fungene_S-Hallin_all_compiled.txt
cat pmoA_hmmer_list3.txt pmoA_hmmer_results_from_nt.txt pmoA_hmmer_results_from_env_nt.txt pmoA_all_fungene.txt pmoa7809.fasta  > pmoA_ncbi_nt_env_nt_fungene_Liebner_all_compiled.txt
cat mmoX_hmmer_list3.txt mmoX_hmmer_results_from_nt.txt mmoX_hmmer_results_from_env_nt.txt mmoX_all_fungene.txt > mmoX_ncbi_nt_env_nt_fungene_all_compiled.txt
cat mcrA_hmmer_list3.txt mcrA_hmmer_results_from_nt.txt mcrA_hmmer_results_from_env_nt.txt mcrA_all_fungene.txt > mcrA_ncbi_nt_env_nt_fungene_all_compiled.txt


#compile all the dataset together:

cat *_all_compiled.txt > all_genes_compiled_ncbi_nt_envnt_fungene.fasta

#cluster sequences to unique ones:

cd-hit-est -i all_genes_compiled_ncbi_nt_envnt_fungene.fasta -o  all_genes_compiled_ncbi_nt_envnt_fungene_unique.fasta -c 1.00

#

# 7°. Submitting Target gene database to MetCap probe design:
#  MetCap probe design:
# all_genes_compiled_ncbi_nt_envnt_fungene_unique.fasta , were submitted to 
# MetCap website (Not available anymore, discuss with Lokesh Manoharan email:lokeshwaran.manoharan@med.lu.se or Dag Ahren email:dag.ahren@med.lu.se in Medical University of Lund ), 
# with 80% clustering and 6 probes per cluster. 

# The database were looked over that no 16S rRNA genes were included into the database (Lokesh? How this was done?)

# 8°. Submitting probes to SeqCap probe order:
# SeqCap probe order:
# produced probe file: Final_N_probes.list1.fasta submitetd to Roche SeqCap probe production:  
# Hyperdesign website to produce probes with Roche SeqCap system 
# https://www.hyperdesign.com/#/
# Probes submitted as Custom probes. 


###############################################################################
# Statitical test for Origical and Targeted reads (Mock community comparison) , plotting the Figure (Script_Fig_S2.sh / Script_ProbeCapture_Fig_S3.sh)
# Statitical test for Shotgun and Targeted reads (Agricultual and wetland soil), plotting the Figure  (Script_Fig_2.sh / Script_Fig_3.sh).



