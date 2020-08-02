
datasets_derived/annotated_sequences/annotated_seq%.csv: datasets_primary/sequencing/hittables/run%_hittable.csv
	@echo Running annotate_sequences.r
	@rscript scripts/annotate_sequences.r $@ $<