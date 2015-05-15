##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Try to determine, based on various parameters, which linker to use for MPI
# applications.
#
ifndef MPI_LD

  $(info ** No MPI linker specified in $$MPI_LD)

  ifdef PE_ENV
    MPI_LD = ftn
  else ifeq '$(OS)' 'AIX'
    MPI_LD = xlf
  else
    MPI_LD = mpif90
  endif

  $(info ** Chose: $(MPI_LD))

endif
