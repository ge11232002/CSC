include $(JKSRC)/inc/common.mk

L = -lm
MYLIBDIR = $(JKSRC)/lib/$(MACHTYPE)
MYLIBS =  $(MYLIBDIR)/jkhgap.a $(MYLIBDIR)/jkweb.a
HG_INC = -I$(JKSRC)/inc -I$(JKSRC)/hg/inc

O = ceScan.o

ceScan: $O $(MYLIBS)
	${CC} ${COPT} ${CFLAGS} -o $(HOME)/bin/ceScan $O $(MYLIBS) $L
	${STRIP} $(HOME)/bin/ceScan${EXE}
