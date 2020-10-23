
annotated_seq%: datasets_primary/sequencing/hittables/run%_hittable.csv
	@echo Running annotate_sequences.r
	@rscript scripts/annotate_sequences.r $@ $<

reference_db: datasets_primary/sequencing/arachnida_reference_seqs.fas datasets_primary/sequencing/odonata_reference_seqs.fas 
	@echo Running mite_reference_db.r
	@rscript scripts/mite_reference_db.r