##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

export ROOT = ../..

include $(ROOT)/make/include.mk
include $(ROOT)/make/compiler.mk
include $(ROOT)/make/fortran/$(FORTRAN_COMPILER).mk

EXTERNAL_LIBRARIES = $(patsubst %,-l%,esmf netcdff netcdf)

FFLAGS += $(OPENMP_ARG)
LDFLAGS += $(OPENMP_ARG)

OBJ_DIR = $(DYNAMO_BUILD_ROOT)/functional-test

ALL_TESTS = $(shell find . -name testframework -prune -o -name "*.py" -print)

SUBJECT_DIR     = $(DYNAMO_BUILD_ROOT)/dynamo
SUBJECT_SUBDIRS = $(shell find $(SUBJECT_DIR) -type d)
SUBJECT_MODINC  = $(patsubst %,-I%,$(SUBJECT_SUBDIRS))

ifdef CRAY_ENVIRONMENT
    STATIC_LIBRARIES = $(EXTERNAL_LIBRARIES)
else
    DYNAMIC_LIBRARIES = $(EXTERNAL_LIBRARIES)
endif

.PHONY: test
test: $(ALL_TESTS)

./%.py: $(OBJ_DIR)/%.x ALWAYS | $(OBJ_DIR)
	@echo -e $(VT_BOLD)Running$(VT_RESET) $@
	$(Q)mpiexec -n 1 ./$@ $<

.PRECIOUS: $(OBJ_DIR)/%.x
$(OBJ_DIR)/%.x: $(OBJ_DIR)/%.o | $(OBJ_DIR)
	@echo -e $(VT_BOLD)Linking$(VT_RESET) $@
	$(Q)$(LDMPI) $(LDFLAGS) $(DYNAMIC_LIBRARIES) -o $@ \
	             $< $(SUBJECT_DIR)/modules.a $(STATIC_LIBRARIES)

.PRECIOUS: $(OBJ_DIR)/%.o
$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)
	@echo -e $(VT_BOLD)Compiling$(VT_RESET) $@
	$(Q)$(FC) $(FFLAGS) -c $(SUBJECT_MODINC) -o $@ $<

$(OBJ_DIR)/%.o: %.F90 | $(OBJ_DIR)
	@echo -e $(VT_BOLD)Compiling$(VT_RESET) $@
	$(Q)$(FC) $(FFLAGS) -c $(SUBJECT_MODINC) -o $@ $<

$(OBJ_DIR):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

.PHONY: ALWAYS
ALWAYS:

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
