include $(JKSRC)/inc/common.mk

L = -lm
MYLIBDIR = $(JKSRC)/lib/$(MACHTYPE)
MYLIBS =  $(MYLIBDIR)/jkhgap.a $(MYLIBDIR)/jkweb.a
HG_INC = -I$(JKSRC)/inc -I$(JKSRC)/hg/inc 

O = chainSubset.o

chainSubset: $O $(MYLIBS)
	${CC} ${COPT} ${CFLAGS} -o $(HOME)/bin/chainSubset $O $(MYLIBS) $L
	${STRIP} $(HOME)/bin/chainSubset${EXE}
