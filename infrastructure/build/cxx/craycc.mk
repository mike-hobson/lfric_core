##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

$(info ** Chosen Cray C++ compiler)

#Older versions of Cray CC compiler don't have --version option
#We can extend this if we need the number in the future like the fortran one
CRAYCC_VERSION_VALID := $(shell CC --version 2> /dev/null || echo "NULLSTRING")

$(info ** Chosen Cray CC compiler version $(CRAYCC_VERSION_VALID))

#Version 015000000 throws errors with cray-c++-rts
ifeq ($(CRAYCC_VERSION_VALID),NULLSTRING)
  CXX_RUNTIME_LIBRARY=cray-c++-rts
else
  CXX_RUNTIME_LIBRARY=
endif