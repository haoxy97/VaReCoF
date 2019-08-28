# ******************* FORTRAN Compiler ********************
FC = g77
FFLAGS = -O3 -pedantic


# ****************** C Compiler *********************
CC = gcc
CFLAGS = -O3 -pedantic

FOR_SOURCES = ch3_h_vdz_corr.f ch_xxx_corr.f cc_000_xxx_corr.f c2h3_ch3_special_corr.f \
              c2h5_h_special_corr.f c2h3_h_special_corr.f pot_aux.f cc_caspt2_corr.f \
              ch_caspt2_corr.f

OBJECTS = $(FOR_SOURCES:.f=.o)

.PHONY: clean tags 

libcorrpot.so : ${OBJECTS} 
	${FC} -shared -o $@ ${OBJECTS}  

clean:
	rm -f *.o
	rm -f *~
	rm -f *.d
	rm -f a.out

tags:
	etags *.cc *.hh *.f

