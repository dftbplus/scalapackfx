############################################################################
#
# Makefile to demonstrate, how to incorporate the library makefile from
# an external makefile by passing the appropriate variables.
#
# Edit "../make.arch" to adapt it to your system.
#
############################################################################

include ../make.arch

# Directories for the necessary libraries
SRCDIR = ../src

############################################################################
# Building the test programs.
#
# You can replace this part with your projects makefile. Make sure, that
# you introduce at least one dependency on the library file (see below).
############################################################################

.SUFFIXES:
.SUFFIXES: .f90 .o


TARGETS = test_desc test_diag test_psyr_pher test_subgrids test_tran

all: $(TARGETS)

# Create dependencies (make sure every targets .o file additionally depends on
# the external libary, as the library modfiles must been created first)
include Makefile.dep
$(TARGETS:=.o) test_common.o: _extlib_scalapackfx

# Include linking rules for targets
define link-target
$(LN) $(LNOPT) -o $@ $(filter-out _%,$^) -L$(SRCDIR) -lscalapackfx $(LIB_SCALAPACK)
endef
include Makefile.targets

%.o: %.f90
	$(FXX) $(FXXOPT) -I$(SRCDIR) -c $<

.PHONY: clean realclean
clean:
	$(MAKE) SRCDIR=. -C $(SRCDIR) -f Makefile.lib clean
	rm -f *.mod *.o _*

distclean: clean
	$(MAKE) SRCDIR=. -C $(SRCDIR) -f Makefile.lib distclean
	rm -f $(TARGETS)

############################################################################
# Invoking the makefile of the librariy to build it in its directory
############################################################################
.PHONY: _FORCED_SUBMAKE_

_extlib_scalapackfx: _FORCED_SUBMAKE_
	touch -r $(SRCDIR)/libscalapackfx.a $@

_FORCED_SUBMAKE_:
	$(MAKE) \
            FXX="$(FXX)" FXXOPT="$(FXXOPT)" \
            LN="$(LN)" LNOPT="$(LNOPT)" \
            M4="$(M4)" M4OPT="$(M4OPT)" \
            SRCDIR=. \
            -C $(SRCDIR) -f Makefile.lib
