.PHONY: clean all

FC = gfortran

OBJS = function_CCSM.o AtmosCCSM.o

target = AtmosCCSM.exe

$(target): $(OBJS)
	$(FC) $(OBJS) liblapack.a librefblas.a -o $@

%.o: %.f90
	$(FC) -c $< -o $@
    
clean:
	-rm -rf *.o *.mod AtmosCCSM.exe tl.bin
