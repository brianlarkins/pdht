#
# main makefile for portals distributed hash table
# 
# author: d. brian larkins
# created: 3/20/2016
#

PORTALS_INCLUDEDIR = /opt/hpctools/include
PORTALS_LIBDIR = /opt/hpctools/lib
PDHT_TOP = .

include pdht.mk

all: pdhtlibs
