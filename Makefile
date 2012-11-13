include ./makefile.include

SUBPROGRAMS=findbasis,findhamiltonian,propagate

$(SUBPROGRAMS):
	$(MAKE) -C $(SUBPROGRAMS)

all:
	$(MAKE) -C $(SUBPROGRAMS) all

.PHONEY: clean
clean::
	$(MAKE) -C $(SUBPROGRAMS) clean
