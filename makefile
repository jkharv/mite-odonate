
annotated_seq%: datasets_primary/sequencing/hittables/run%_hittable.csv
	@echo Running annotate_sequences.r
	@rscript scripts/annotate_sequences.r $@ $<