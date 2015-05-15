##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Set up a few bits and bobs relating to the build system.
# This file is intended to be "include"d from other make files.
# It expects ROOT to point to the root of the project tree. i.e. the directory
# containing the top level Makefile.
#
BUILD_DIR = $(ROOT)/build
TOOL_DIR = $(ROOT)/tools
MAKE_DIR = $(ROOT)/make

OS = $(shell uname -s)
MACHINE = $(shell uname -m)

ifdef VERBOSE
    Q :=
else
    Q = @
endif

export Q

VT_BOLD = \\x1b[1m
VT_RESET = \\x1b[0m
