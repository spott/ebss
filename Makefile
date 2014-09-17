include ./makefile.include

SUBPROGRAMS=findbasis findhamiltonian propagate

findbasis:
	$(MAKE) -C findbasis/

findhamiltonian:
	$(MAKE) -C $@

propagate:
	$(MAKE) -C $@

all: $(SUBPROGRAMS)

.PHONY: all
