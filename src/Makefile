.SUFFIXES:
.SUFFIXES: .f90 .o

RM     = /bin/rm -f

FC     = gfortran -ffree-form -fcheck=all -Wunused -Wuninitialized
#FC     = gfortran -ffree-form -fopenmp
FC_old = gfortran -std=legacy -ffixed-line-length-none

%.o: %.f90
	$(FC) $< -c

%.o: %.f
	$(FC_old) $< -c

$(DYNBAS)/bin/dynbas: dynbas.f90 dynbas_pieces_module.o nnls.o brent.o locate.o
	$(FC) dynbas.f90 dynbas_pieces_module.o nnls.o  brent.o locate.o -o dynbas
	if [ ! -d "$(DYNBAS)/bin" ]; then mkdir $(DYNBAS)/bin; fi
	mv dynbas $(DYNBAS)/bin/.

brent.o: brent.f
locate.o: locate.f
nnls.o: nnls.f
dynbas_pieces_module.o: dynbas_pieces_module.f90 nnls.o

install: $(DYNBAS)/bin/dynbas
	make uninstall
	if [ -d "$(HOME)/.bin" ]; then ln -s $(DYNBAS)/bin/dynbas $(HOME)/.bin/dynbas; fi
	if [ -d "$(HOST_HOME)/.bin" ]; then ln -s $(DYNBAS)/bin/dynbas $(HOST_HOME)/.bin/dynbas; fi

uninstall:
	if [ -e "$(HOME)/.bin/dynbas" ]; then $(RM) $(HOME)/.bin/dynbas; fi
	if [ -e "$(HOST_HOME)/.bin/dynbas" ]; then $(RM) $(HOST_HOME)/.bin/dynbas; fi

clean:
	$(RM) *.o *.mod
	if [ -d "$(DYNBAS)/bin" ]; then $(RM) $(DYNBAS)/bin/*; fi
	make uninstall
