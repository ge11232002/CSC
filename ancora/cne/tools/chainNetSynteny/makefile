include $(JKSRC)/inc/common.mk

L = -lm $(MYSQLLIBS)
MYLIBDIR = $(JKSRC)/lib/$(MACHTYPE)
MYLIBS =  $(MYLIBDIR)/jkhgap.a $(MYLIBDIR)/jkweb.a
HG_INC = -I$(JKSRC)/inc -I$(JKSRC)/hg/inc

O = chainNetSynteny.o

chainNetSynteny: $O $(MYLIBS)
	${CC} ${COPT} ${CFLAGS} -o $(HOME)/bin/chainNetSynteny $O $(MYLIBS) $L
	${STRIP} $(HOME)/bin/chainNetSynteny${EXE}
