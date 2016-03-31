degeneracy_tables.R:
	python makeprimertables.py

loaddata: degeneracy_tables.R
	R CMD BATCH makeprimertables.R

all: loaddata

clean:
	rm degeneracy_tables.R .RData

