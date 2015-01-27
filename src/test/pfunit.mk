##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

CMAKE ?= cmake

PFUNIT_SOURCE_DIR = $(abspath ../pfunit)
PFUNIT_BUILD_DIR = $(BUILD_DIR)/pfunit

include $(ROOT)/include.mk

COMPILER_NAME = $(shell basename $(FC))
ifeq '$(COMPILER_NAME)' 'ifort'
  PFUNIT_COMPILER_ID = Intel
  ifeq ($(shell test $(IFORT_VERSION) -lt 0130000; echo $$?), 0)
    $(error pFUnit will only compile with ifort v13 or later.)
  else ifeq ($(shell test $(IFORT_VERSION) -lt 0140000; echo $$?), 0)
    export CPPFLAGS += -DINTEL_13
    export FPPFLAGS += -DINTEL_13
  endif
  FFLAGS := $(filter-out -warn errors all,$(FFLAGS)) # pFUnit isn't that well put together
else ifeq '$(COMPILER_NAME)' 'gfortran'
  PFUNIT_COMPILER_ID = GNU
  ifeq ($(shell test $(GFORTRAN_VERSION) -lt 040500), 0)
    $(error pFUnit will only compile with gfortran v4.5 or later.)
  endif
  FFLAGS := $(filter-out -Werror,$(FFLAGS)) # pFUnit isn't that well put together
else ifeq '$(COMPILER_NAME)' 'nagfor'
  PFUNIT_COMPILER_ID = NAG
else ifeq '$(COMPILER_NAME)' 'xlf'
  PFUNIT_COMPILER_ID = XL
else
  $(error Unrecognised compiler "$(FC)")
endif

DRIVER_DIR = $(dir $(DRIVER_OBJ))

$(DRIVER_OBJ): $(PFUNIT_INSTALL_DIR)/include/driver.F90 $(DRIVER_DIR)testSuites.inc | $(DRIVER_DIR)
	@echo Compiling $@
	$(Q)$(FC) $(FFLAGS) -c -I$(PFUNIT_INSTALL_DIR)/mod -I$(DRIVER_DIR) \
	          -DBUILD_ROBUST -o $@ $<

$(PFUNIT_INSTALL_DIR)/include/driver.F90: $(PFUNIT_BUILD_DIR)/Makefile
	@echo Building pFUnit
	$(Q)$(MAKE) -C $(PFUNIT_BUILD_DIR)
	$(Q)$(MAKE) -C $(PFUNIT_BUILD_DIR) tests install

$(PFUNIT_BUILD_DIR)/Makefile: | $(PFUNIT_BUILD_DIR)
	@echo Configuring pFUnit
	$(Q)cd $(PFUNIT_BUILD_DIR); $(CMAKE) -DINSTALL_PATH=$(PFUNIT_INSTALL_DIR) \
	                                     $(PFUNIT_SOURCE_DIR)

$(PFUNIT_BUILD_DIR) $(dir $(DRIVER_OBJ)):
	@echo Creating $@
	$(Q)mkdir $@

ALL_TESTS = $(shell find . -name "*.pf")

$(DRIVER_DIR)testSuites.inc: $(DRIVER_DIR)testSuites.inc.new $(ALL_TESTS)
	@echo Replacing $@ with $<
	@mv -f $< $@

$(DRIVER_DIR)testSuites.inc.new: ALWAYS | $(DRIVER_DIR)
	@echo Clearing $@
	@echo > $@

$(ALL_TESTS): ALWAYS
	@echo Adding $@
	@echo ADD_TEST_SUITE\($(notdir $(basename $@))_suite\) >> $(DRIVER_DIR)testSuites.inc.new

.PHONY: ALWAYS
ALWAYS:

.PHONY: clean
clean:
	-rm -rf $(PFUNIT_BUILD_DIR)
	-rm -rf $(PFUNIT_INSTALL_DIR)
