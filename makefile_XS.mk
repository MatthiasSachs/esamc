CC = gcc # compiler used for compilation of c-code
CPPC = g++ # compiler used for compilation of c++ code
GSL_LIB_PATH = /home/xshang/mypref/lib/ # Path to directory where libgsl.a is located
HDF5_LIB_PATH = /home/xshang/mypref/lib/ # Path to directory where libhdf5_hl.a and libhdf5.a are located
GSL_HEADER_PATH = /home/xshang/mypref/include/ # Path to directory where the directory "gsl/" containing all the gsl *.h files is located
HDF5_HEADER_PATH = /home/xshang/mypref/include/ # Path to directory where all H5*.h files are located
CFLAGS_OPTIONAL = -std=c99 #optional flags for c compiler
CPPFLAGS_OPTIONAL = -std=c++11 -o esamc # optional flags for c++ compiler

SRCS := $(wildcard *.cpp)
CFLAGS_DEFAULT = -I$(GSL_HEADER_PATH) -I$(HDF5_HEADER_PATH) -L$(GSL_LIB_PATH) -L$(GSL_LIB_PATH) -lgsl -lgslcblas -lhdf5 -lhdf5_hl
CPPFLAGS_DEFAULT = -I$(GSL_HEADER_PATH) -I$(HDF5_HEADER_PATH) -L$(GSL_LIB_PATH) -L$(GSL_LIB_PATH) -lgsl -lgslcblas -lhdf5 -lhdf5_hl   
all: precompile_c_files compile_main_code

precompile_c_files: 
	${CC} expokit_translation2.c -c ${CFLAGS_OPTIONAL} ${CFLAGS_DEFAULT}
compile_main_code:
	${CPPC} $(SRCS) ${CPPFLAGS_OPTIONAL} ${CFLAGS_DEFAULT} ./expokit_translation2.o
clean:
	@echo "Cleaning up..."
	@rm *.o

#
# $ make -f makefile_XS.mk 
# $ LD_LIBRARY_PATH=/home/xshang/mypref/lib
# $ export LD_LIBRARY_PATH
# $ ./esamc
#
# https://www.gnu.org/software/gsl/doc/html/usage.html#f4
#
# https://stackoverflow.com/questions/4470838/g-linking-issue-with-gsl
#