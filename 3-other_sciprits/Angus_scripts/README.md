# README

This readme provides a brief overview of the function of the different scripts used
in the probe capture project.

scripts
    |- add_ko_col.py        -   Adds the specified KO (arg 2) identifier to the provided TSV (arg 1) as the last column.
    |- best_ko_mappings.py  -   Gets the KO value for each protein in the provided TSV with the best e-value (assumes e-values are in the 5th column)
    |- best_non_overlapping_hits.py
                            -   Gets the top two non-overlapping hits (where overlap is defined as no more than 5 amino acids) per query in a BLAST file produced
                                by Diamond BLAST (Buchfink et al. 2021)
    |- build_best_read_mappings.py
                            -   Adds KO identifiers from full length proteins to which reads mapped, used to determine trye values of read-mappings. 
    |- build_cult_fwd_lookup.py
                            -   Builds a lookup table to map reads to their respective proteins (based on genomes in the culture) using a blast table (i.e., 
                                reads against the genomes.)
    |- build_results.py     -   Builds the results table, mapping reads to their respective proteins via the read-to-protein mapping file, used to check if
                                the KO's match in downstream analyses.
    |- check_overlaps.py    -   Checks for potential gene fusions by checking for overlapping results in the domains file.
    |- check_resulst.py     -   Validation script, used to count matches to accessions to inspect results for errors (e.g., in case a results has more than
                                4 hits to it)
    |- count_incorrect_kos.py
                            -   Checks if KOs do not match in the results file to identify incorrectly called reads.
    |- count_stats.py       -   Generates statistics for each provided input file (arg 3...), each representing the results for a given gene. Results are based
                                on true positive mapping (read to its associated KO; arg 1) and counts for all positive results (arg 2). Results include the
                                true positive, true negatives, false negative, precision, and recall.
    |- dedup.py             -   Short script to remove duplicate reads from a FASTA file.
    |- generate_mappings.py -   Reads a BLAST results table and gets the top four hits per read, accounting for all four possible reading frames (up to 2 hits
                                forward and 2 hits in reverse.)
    |- get_best_hits.py     -   Gets the best HMM hits per read. This includes up to two non-overlapping hits for each protein that was checked (accounts for
                                possible gene fusions)
    |-  get_best_map_per_read.py
                            -   Gets the best result for every read in a table (based on minimum score in the 11th column)
    |-  get_missing.py      -   Gets missing results by getting all identifiers in arg 1 and checking if they exist in arg 2 (identifiers from first column)
    |-  join_kos.py         -   Creates a lookup of reads to KOs by checking a mapping file. Results will be in a TSV with the first column representing the reads
                                and the second column being a comma separated list of all KOs found in the input table (arg 1) matching to that KO.
    |-  map_final_results.py
                            -   Generates a mapping of reads to KOs via protein mappings. The inputs are, e.g., the top 2 BLAST results, HMM results for the
                                genomes, and the actual results/best HMM per frame. 
    |-  map_final_resultsDOM.py
                            -   Same as above, includes domain information in the output.
    |-  match_results.py    -   Used to check the results from the in-house HMM files (provided by Henri Siljanen) against the predicted true KOs (from protein 
                                mappings.)
    |-  build_project.sh    -   Used to build main project directory and generate basic files for server-side analysis (usded in the run*.sh scripts for slurm)
    |-  run_exec.sh         -   Runs exec_annotation, used to get the KOFAM results (KofamScan: https://github.com/takaram/kofam_scan) using a slurm batch job
    |-  run_hmm.sh          -   Runs hmmseach using a slurm batch job
    |-  run_kofam.sh        -   Runs exec_annotation on idividual *.faa files, primarily used for the reads. Uses a slurm batch job.
    |- filters/             -   This directory contains different scripts that were used to filter data tables based on various parameters (see files for further 
                                information)