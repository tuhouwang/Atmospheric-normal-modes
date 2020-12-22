.PHONY: clean all

FC = gfortran

OBJS = function_CTSM.o AtmosCTSM.o

target = AtmosCTSM.exe

$(target): $(OBJS)
	$(FC) $(OBJS) liblapack.a librefblas.a -o $@

%.o: %.f90
	$(FC) -c $< -o $@
    
clean:
	-rm -rf *.o *.mod AtmosCTSM.exe tl.bin
