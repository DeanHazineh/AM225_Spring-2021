# Include the configuration file from the parent directory
include ../config.mk

# Lists of files to be built
objs = StrassenFunctions.o  conj_grad.o file_output.o rbf.o lp_solve.o rbfAnalysis.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs = Strassen_Driver rbf_test rbf_test2 rbf_time rbf_driver

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) $(fftw_iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

Strassen_Driver: Strassen_Driver.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

rbf_driver: rbf_driver.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

.PHONY: clean depend
