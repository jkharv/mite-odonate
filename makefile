#seqtab_otu.csv seqtab_asv.csv
assign_taxonomy:
	@echo Running assign_taxonomy.r
	@Rscript scripts/assign_taxonomy.r seqtab_asv.csv
	
network:
	@Rscript scripts/make_network.r
