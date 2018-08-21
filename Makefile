#
# main makefile for portals distributed hash table
# 
# author: d. brian larkins
# created: 3/20/2016
#

PORTALS_INCLUDEDIR = $(HOME)/opt/include
PORTALS_LIBDIR = $(HOME)/opt/lib
PDHT_TOP = .

include pdht.mk

all: pdhtlibs
