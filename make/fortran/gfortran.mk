##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Determine the version of the Intel Fortran compiler.
#
# This macro is evaluated now (:= syntax) so it may be used as many times as
# desired without wasting time rerunning it.
#
GFORTRAN_VERSION     := $(shell $(FC) -dumpversion 2>&1 \
                        | awk -F . '{ printf "%02i%02i%02i", $$1, $$2, $$3 }')
$(info ** GNU Fortran version $(GFORTRAN_VERSION))

F_MOD_DESTINATION_ARG     = -J
