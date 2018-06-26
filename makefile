F90 = \ifort 
F90omp = \ifort 
a.out : global_constants.o Woods_Saxon.o
	${F90} global_constants.o Woods_Saxon.o
global_constants.o : global_constants.f90
	${F90omp} -c global_constants.f90
Woods_Saxon.o : global_constants.f90 Woods_Saxon.f90
	${F90omp} -c Woods_Saxon.f90

.PHONY : clean
clean :
	$(RM) *.o *.mod a.out gnu EPS_Files/* Waves/No_LS/* Waves/With_LS/*
	$(RM) TEX_Files/Woods_Saxon.dvi  TEX_Files/Woods_Saxon.pdf
