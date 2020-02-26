CC = g++
LD = g++

#EIGEN_INCLUDE = -I/home/ttaka/eigen-3.3.4
EIGEN_INCLUDE = -I/home/ttaka/eigen-3.3.7
REDSVD_INCLUDE = -I.

# make sure the corresponding modules have been added
# and use 'extern "C"{} for declaring lapack routines'
# LDPATH = -L $(BLAS_DIR) -L $(LAPACK_DIR) -L $(FFTW_DIR)
#LDPATH = -L/usr/
LDPATH =

#LDFLAGS = -llapack -lblas

MKLROOT = /opt/intel/mkl
LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LDFLAGS += -Wl,-rpath,${MKLROOT}/lib/intel64

LDFLAGS += -lgfortran
LDFLAGS += -lgomp -lpthread
LDFLAGS += -lm


#ttaka# FLAGS  = -g -msse2 -Wall -O2 -I $(BLAS_INCLUDE) -I /Users/ericdarve/Documents/Pieter/eigenorig
FLAGS = -Wall
FLAGS += -O2 -mtune=native -march=native $(BLAS_INCLUDE) ${REDSVD_INCLUDE} ${EIGEN_INCLUDE}
FLAGS += -g -D_DEBUG -DMYDEBUG

PFLAG  =

FLAGS += -DEIGEN_USE_MKL_ALL
#FLAGS += -DEIGEN_DONT_PARALLELIZE

FLAGS += -fopenmp

#FLAGS += -DPARA_CHECK_RANK
#FLAGS += -DIFMM_PARALLELIZE

#FLAGS += -DDISABLE_JACOBISVD


PROG = bbifmm
OBJS = main.o ifmm.o set_L2P_operator.o set_L2L_operator.o set_M2M_operator.o set_P2M_operator.o set_M2P_operator.o set_P2L_operator.o set_M2L_operator.o set_P2P_operator.o set_nDOF.o set_RHS.o set_xyz.o set_points.o check_rank.o get_weights.o transfer_M2L_to_P2P.o transfer_M2M_to_P2M.o transfer_L2L_to_L2P.o compute_Tk.o transfer_RHS.o ACA_FullyPivoted.o
OBJS += parallel.o
OBJS += timer.o

all:    $(PROG)

main.o: main.cpp ifmm.h
	$(CC) $(FLAGS) -c main.cpp $(PFLAG)

ifmm.o: ifmm.cpp ifmm.h container.h
	$(CC) $(FLAGS) -c ifmm.cpp $(PFLAG)

set_L2P_operator.o: set_L2P_operator.cpp set_L2P_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_L2P_operator.cpp $(PFLAG)

set_L2L_operator.o: set_L2L_operator.cpp set_L2L_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_L2L_operator.cpp $(PFLAG)

set_M2M_operator.o: set_M2M_operator.cpp set_M2M_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_M2M_operator.cpp $(PFLAG)

set_P2M_operator.o: set_P2M_operator.cpp set_P2M_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_P2M_operator.cpp $(PFLAG)

set_M2P_operator.o: set_M2P_operator.cpp set_M2P_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_M2P_operator.cpp $(PFLAG)

set_P2L_operator.o: set_P2L_operator.cpp set_P2L_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_P2L_operator.cpp $(PFLAG)

set_M2L_operator.o: set_M2L_operator.cpp set_M2L_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_M2L_operator.cpp $(PFLAG)

set_P2P_operator.o: set_P2P_operator.cpp set_P2P_operator.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_P2P_operator.cpp $(PFLAG)

set_nDOF.o: set_nDOF.cpp set_nDOF.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_nDOF.cpp $(PFLAG)

set_RHS.o: set_RHS.cpp set_RHS.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_RHS.cpp $(PFLAG)

set_xyz.o: set_xyz.cpp set_xyz.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_xyz.cpp $(PFLAG)

set_points.o: set_points.cpp set_points.h container.h ifmm.h
	$(CC) $(FLAGS) -c set_points.cpp $(PFLAG)

check_rank.o: check_rank.cpp check_rank.h container.h ifmm.h
	$(CC) $(FLAGS) -c check_rank.cpp $(PFLAG)

get_weights.o: get_weights.cpp get_weights.h container.h ifmm.h
	$(CC) $(FLAGS) -c get_weights.cpp $(PFLAG)

compute_Tk.o: compute_Tk.cpp compute_Tk.h container.h
	$(CC) $(FLAGS) -c compute_Tk.cpp $(PFLAG)

transfer_M2L_to_P2P.o: transfer_M2L_to_P2P.cpp transfer_M2L_to_P2P.h container.h ifmm.h
	$(CC) $(FLAGS) -c transfer_M2L_to_P2P.cpp $(PFLAG)

transfer_M2M_to_P2M.o: transfer_M2M_to_P2M.cpp transfer_M2M_to_P2M.h container.h ifmm.h
	$(CC) $(FLAGS) -c transfer_M2M_to_P2M.cpp $(PFLAG)

transfer_L2L_to_L2P.o: transfer_L2L_to_L2P.cpp transfer_L2L_to_L2P.h container.h ifmm.h
	$(CC) $(FLAGS) -c transfer_L2L_to_L2P.cpp $(PFLAG)

transfer_RHS.o: transfer_RHS.cpp transfer_RHS.h container.h ifmm.h
	$(CC) $(FLAGS) -c transfer_RHS.cpp $(PFLAG)

ACA_FullyPivoted.o: ACA_FullyPivoted.cpp ACA_FullyPivoted.h container.h ifmm.h
	$(CC) $(FLAGS) -c ACA_FullyPivoted.cpp $(PFLAG)

parallel.o: parallel.cpp container.h ifmm.h
	$(CC) $(FLAGS) -c parallel.cpp $(PFLAG)

$(PROG): $(OBJS)
	$(LD) -o $(PROG) $(OBJS) $(LDPATH) $(LDFLAGS) $(PFLAG)

clean:
	/bin/rm -rf *.o *~ \#*\# $(PROG) *core

tar:
#	tar cvfz 3DbbiFMM.tar.gz Makefile main.cpp ifmm.cpp ifmm.h container.h
	tar cvfz tmp.tar.gz Makefile *.cpp *.h *.hpp
