##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Compiler specific details: NAG Fortran and GCC C/C++
#
$(info ** Chosen NAG Fortran compiler)

FFLAGS_COMPILER           =
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -O4
FFLAGS_DEBUG              = -g -gline -colour
FFLAGS_INIT               = -nan
FFLAGS_RUNTIME            = -C=all -C=undefined -mtrace=all -nan

LDFLAGS_COMPILER =

