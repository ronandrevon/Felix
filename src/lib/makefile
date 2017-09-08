#/usr/bin/bash

include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk

DIRCIFTBX=../ciftbx/
DIRQUADPACK=../quadpack/
DIRFELIX=../felix/

LIBFELIX=$(MYPLATFORM)/lib$(PRECISION)felix.a
LIBQUADPACK=$(MYPLATFORM)/lib$(PRECISION)quadpack.a
LIBCIFTBX=$(MYPLATFORM)/lib$(PRECISION)ciftbx.a

ourlibs: $(LIBFELIX) $(LIBQUADPACK) $(LIBCIFTBX)


$(LIBFELIX): $(FELIX) 
	mkdir -p $(MYPLATFORM)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(LIBQUADPACK):  $(QUADPACK)
	mkdir -p $(MYPLATFORM)
	$(AR) $(ARFLAGS) $@ $<
	$(RL) $(RLFLAGS) $@

$(LIBCIFTBX):  $(CIFTBX)
	mkdir -p $(MYPLATFORM)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

