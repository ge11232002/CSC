include $(JKSRC)/inc/common.mk	

L = -lm
MYLIBDIR = $(JKSRC)/lib/$(MACHTYPE)
MYLIBS =  $(MYLIBDIR)/jkhgap.a $(MYLIBDIR)/jkweb.a
HG_INC = -I$(JKSRC)/inc -I$(JKSRC)/hg/inc 

O = netChainSubsetRelaxed.o

netChainSubsetRelaxed: $O $(MYLIBS)
	${CC} ${COPT} ${CFLAGS} -o $(HOME)/bin/netChainSubsetRelaxed $O $(MYLIBS) $L
	${STRIP} $(HOME)/bin/netChainSubsetRelaxed${EXE}
