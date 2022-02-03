# Makefile for the DLPFC-SCZ Project

RS=Rscript
DATA_DIR=data/
RES_DIR=results/
FIG_DIR=figures/
MAP_DIR=maps/

.PHONY: complete

complete: clean all

clean:
	rm -rfv $(DATA_DIR)
	rm -rfv $(RES_DIR)
	rm -rfv $(FIG_DIR)
	rm -rfv $(MAP_DIR)
	$(RS) 0-setup.R

make_lists:

process_drugs:

process_groups:

process_diseases:

analyse:

visualise:

all: make_lists process_drugs process_groups process_diseases analyse visualise
