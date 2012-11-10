include ./makefile.include

SUBPROGRAMS=findbasis,findhamiltonian,propagate

$(subprograms):
	$(MAKE) -C $(subprograms)

#.PHONEY: clean
#clean:
	#$(MAKE) -C $(subprograms) clean
