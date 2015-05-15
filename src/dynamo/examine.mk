##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
#
# Python's support for SQLite seems to hobble its concurrency support so this
# stuff has to be done serially.

.NOTPARALLEL:

ROOT = ../..
include $(ROOT)/make/include.mk

DATABASE = $(OBJ_DIR)/dependencies.db

ALL_SRC = $(shell find . -name "*.[Ff]90")
TOUCH_FILES = $(patsubst ./%.f90,$(OBJ_DIR)/%.t,$(patsubst ./%.F90,$(OBJ_DIR)/%.t,$(ALL_SRC)))
IGNORE_ARGUMENTS := $(patsubst %,-ignore %,$(IGNORE_DEPENDENCIES))

$(OBJ_DIR)/programs.mk: $(OBJ_DIR)/dependencies.mk | $(OBJ_DIR)
	@echo -e $(VT_BOLD)Building$(VT_RESET) $@
	$(Q)$(TOOL_DIR)/ProgramObjects -database $(DATABASE) $@

$(OBJ_DIR)/dependencies.mk: $(TOUCH_FILES) | $(OBJ_DIR) $(OBJ_SUBDIRS)
	@echo -e $(VT_BOLD)Building$(VT_RESET) $@
	$(Q)$(TOOL_DIR)/DependencyRules -database $(DATABASE) $@

$(OBJ_DIR)/%.t: %.F90 | $(OBJ_DIR) $(OBJ_SUBDIRS)
	@echo -e $(VT_BOLD)Analysing$(VT_RESET) $<
	$(Q)$(TOOL_DIR)/DependencyAnalyser $(IGNORE_ARGUMENTS) \
	                                   $(DATABASE) $< && touch $@

$(OBJ_DIR)/%.t: %.f90 | $(OBJ_DIR) $(OBJ_SUBDIRS)
	@echo -e $(VT_BOLD)Analysing$(VT_RESET) $<
	$(Q)$(TOOL_DIR)/DependencyAnalyser $(IGNORE_ARGUMENTS) \
	                                   $(DATABASE) $< && touch $@

$(OBJ_DIR) $(OBJ_SUBDIRS):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@
