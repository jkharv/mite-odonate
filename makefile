# Author Jake Harvey - jakekharvey@gmail.com

datasets_derived/sequencing/mite_sequences_annotated.csv: datasets_derived/sequencing/mite_sequences.csv scripts/assign_taxonomy.r | datasets_derived/sequencing/reference_db.fasta
	Rscript scripts/assign_taxonomy.r
	echo Assignning taxonomy to sequences.

datasets_derived/bin_network.csv datasets_derived/prob_network.csv: datasets_derived/sequencing/mite_sequences_annotated.csv scripts/make_network.r
	Rscript scripts/make_network.r
	echo Making network.

datasets_derived/mite_phylo_scale.csv: datasets_derived/bin_network.csv scripts/indices.r
	Rscript scripts/indices.r
	echo Calculating specialization indices.

datasets_derived/odonate_summaries.csv: scripts/odonate_summaries.r datasets_derived/bin_network.csv datasets_derived/mite_phylo_scale.csv
	Rscript scripts/odonate_summaries.r
	echo Summarizing data at the odonate species level.

datasets_derived/mite_summaries.csv: scripts/odonate_summaries.r datasets_derived/bin_network.csv datasets_derived/mite_phylo_scale.csv
	Rscript scripts/odonate_summaries.r
	echo Summarizing data at the odonate species level.
