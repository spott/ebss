include ./makefile.include

SUBPROGRAMS=findbasis findhamiltonian propagate

nonlinear:
	$(MAKE) -C utils nonlinear

utils:
	$(MAKE) -C utils all

all: $(SUBPROGRAMS) utils

$(SUBPROGRAMS):
	$(MAKE) -C $@ $@

clean::
	$(foreach prog, $(SUBPROGRAMS), $(MAKE) -C $(prog) clean;)
	$(MAKE) -C utils clean;

.PHONY: all clean $(SUBPROGRAMS) nonlinear utils
