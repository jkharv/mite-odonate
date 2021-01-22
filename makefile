# Author Jake Harvey - jakekharvey@gmail.com

# datasets_primary/
# 	phylogeny/
# 		FinalOdonatetree.tre
#	sequencing/
#		raw_fastq/
#		arachnida_reference_seqs.fas
#		odonata_reference_seqs.fas
#	mite_samples.csv
#
# datasets_derived/
#	sequencing/
#		mite_sequences.csv
#		mite_sequences_annotated.csv
#		mite_sequences_otu.csv
#		reference_db.fasta
#	bin_network.csv
#	prob_network.csv
#	controls.csv
#	mite_phylo_scale.csv
#
# scripts/
#	assign_taxonomy.r
#	clean_sequences.r
#	indices.r
#	install_dependencies.r
#	make_network.r
#	mite_reference_db.r
#	phylo_scale.r
#	figures/
#		phylo_mite_scale.r
#		phylo_scale_histogram.r

datasets_derived/sequencing/mite_sequences.csv datasets_derived/sequencing/mite_sequences_otu.csv: datasets_primary/sequencing/raw_fastq/*.fastq.gz scripts/clean_sequences.r
	Rscript scripts/clean_sequences.r
	echo Cleaning raw sequence data

datasets_derived/sequencing/reference_db.fasta: datasets_primary/sequencing/arachnida_reference_seqs.fas datasets_primary/sequencing/odonata_reference_seqs.fas
	Rscript scripts/mite_reference_db.r
	echo Creating reference database for mite sequences

datasets_derived/sequencing/mite_sequences_annotated.csv: mite_sequences.csv scripts/assign_taxonomy.r | datasets_derived/sequencing/reference_db.fasta
	Rscript scripts/assign_taxonomy
	echo Assignning taxonomy to sequences

datasets_derived/bin_network.csv datasets_derived/prob_network.csv: datasets_derived/sequencing/mite_sequences_annotated.csv scripts/make_network.r
	Rscript scripts/make_network.r
	echo Making network

datasets_derived/mite_phylo_scale.csv: datasets_derived/bin_network.csv scripts/indices.r scripts/phylo_scale.r
	Rscript scripts/indices.r
	echo Calculating phylogenetic scale and specialization indices

figures/phylo_plot.svg: datasets_derived/mite_phylo_scale.csv scripts/figures/phylo_mite_scale.r
	Rscript scripts/figures/phylo_mite_scale.r

figures: figures/phylo_plot.svg