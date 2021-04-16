# Include the configuration file from the parent directory
include ../config.mk

# Lists of files to be built
objs = MRI_DST_FFTW.o MRI_BandedFD.o file_output.o MRI_DST_FFTW_V2.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs = MRI_Driver banded_test

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

MRI_Driver: MRI_Driver.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

banded_test: banded_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)


.PHONY: clean depend
