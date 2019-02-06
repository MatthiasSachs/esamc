CC = gcc # compiler used for compilation of c-code
CPPC = g++ # compiler used for compilation of c++ code
GSL_LIB_PATH = /usr/local/Cellar/gsl/2.4/lib/ # Path to directory where libgsl.a is located
HDF5_LIB_PATH = /usr/local/lib/ # Path to directory where libhdf5_hl.a and libhdf5.a are located
GSL_HEADER_PATH = /usr/local/Cellar/gsl/2.4/include/ # Path to directory where the directory "gsl/" containing all the gsl *.h files is located
HDF5_HEADER_PATH = /usr/local/lib/ # Path to directory where all H5*.h files are located
CFLAGS_OPTIONAL =  #optional flags for c compiler 
CPPFLAGS_OPTIONAL = -o esamc # optional flags for c++ compiler 

SRCS := $(wildcard *.cpp)
CFLAGS_DEFAULT = -I$(GSL_HEADER_PATH) -I$(HDF5_HEADER_PATH) -L$(GSL_LIB_PATH) -L$(GSL_LIB_PATH) -lgsl -lhdf5 -lhdf5_hl
CPPFLAGS_DEFAULT = -I$(GSL_HEADER_PATH) -I$(HDF5_HEADER_PATH) -L$(GSL_LIB_PATH) -L$(GSL_LIB_PATH) -lgsl -lhdf5 -lhdf5_hl   
all: precompile_c_files compile_main_code

precompile_c_files: 
	${CC} expokit_translation2.c -c ${CFLAGS_OPTIONAL} ${CFLAGS_DEFAULT}
compile_main_code:
	${CPPC} $(SRCS) ${CPPFLAGS_OPTIONAL} ${CFLAGS_DEFAULT} ./expokit_translation2.o
clean:
	@echo "Cleaning up..."
	@rm *.o


