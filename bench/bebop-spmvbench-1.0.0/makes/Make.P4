# CITRIS-specific makefile for standalone SpMV benchmark
# author: hormozd

# Basic commands
CD           = cd
RM           = /bin/rm -f

# We're on this platform
ARCH         = $(arch)

# Compiler
CC           = gcc
CFLAGS       = -O3 -Wall -std=c99
LINKER       = gcc
LINKFLAGS    = 
ARCHIVER     = ar
ARFLAGS      = rc
RANLIB       = ranlib
