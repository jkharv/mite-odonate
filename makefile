# Author Jake Harvey - jakekharvey@gmail.com

datasets_derived/sequencing/mite_sequences.csv datasets_derived/sequencing/mite_sequences_otu.csv: datasets_primary/sequencing/raw_fastq/*.fastq.gz scripts/clean_sequences.r
	Rscript scripts/clean_sequences.r
	echo Cleaning raw sequence data

datasets_derived/sequencing/reference_db.fasta: datasets_primary/sequencing/arachnida_reference_seqs.fas datasets_primary/sequencing/odonata_reference_seqs.fas
	Rscript scripts/mite_reference_db.r
	echo Creating reference database for mite sequences

datasets_derived/sequencing/mite_sequences_annotated.csv: datasets_derived/sequencing/mite_sequences.csv scripts/assign_taxonomy.r | datasets_derived/sequencing/reference_db.fasta
	Rscript scripts/assign_taxonomy.r
	echo Assignning taxonomy to sequences

datasets_derived/bin_network.csv datasets_derived/prob_network.csv: datasets_derived/sequencing/mite_sequences_annotated.csv scripts/make_network.r
	Rscript scripts/make_network.r
	echo Making network

datasets_derived/mite_phylo_scale.csv: datasets_derived/bin_network.csv scripts/indices.r scripts/phylo_scale.r
	Rscript scripts/indices.r
	echo Calculating phylogenetic scale and specialization indices

datasets_derived/odonate_summaries.csv: scripts/odonate_summaries.r datasets_primary/2015_data.csv datasets_primary/2019_data.csv datasets_primary/2020_data.csv datasets_derived/bin_network.csv datasets_derived/mite_phylo_scale.csv
	Rscript scripts/odonate_summaries.r
	echo Calculating Odonate abundances and mite prevalences
