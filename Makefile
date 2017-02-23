include ./makefile.include

SUBPROGRAMS=findbasis findhamiltonian propagate

findbasis:
	$(MAKE) -C $@ $@

findhamiltonian:
	$(MAKE) -C $@ $@

propagate:
	$(MAKE) -C $@ $@

all: $(SUBPROGRAMS)

clean::
	$(foreach prog, $(SUBPROGRAMS), $(MAKE) -C $(prog) clean;)

.PHONY: all clean $(SUBPROGRAMS)
