##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Compiler specific details: IBM XL Fortran
#
$(info ** Chosen IBM XL Fortran compiler)

FFLAGS_COMPILER           = -qxlf2003=polymorphic
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -O4
FFLAGS_DEBUG              = -g -qdbgfmt=dwarf
FFLAGS_WARNINGS           = -qinfo=all -qhalt=w
FFLAGS_INIT               = -qinitalloc=ffffffff -qinitauto=ffffffff
FFLAGS_RUNTIME            = -qddim -qstackprotect=all -qfloat=nans -qflttrap=enable:invalid:nanq:overflow:underflow:zerodivide
