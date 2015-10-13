OS=$(shell lsb_release -si)
ARCH=$(shell uname -m | sed 's/x86_//;s/i[3-6]86/32/')
VER=$(shell lsb_release -sr)

CFLAGS=-c -std=c++0x -O3 -pthread -flto -DNDEBUG -Wall -Wpointer-arith -Wextra -Wno-literal-suffix #-Wshadow -pedantic -Wcast-qual -Wstrict-overflow=5 -m64 -Ofast  -march=native -funroll-loops
LDFLAGS= -pthread -flto -lm -lgsl -lgslcblas -lsqlite3

ifeq ($(strip $(OS)),$(filter $(strip $(OS)),SUSE LINUX openSUSE project)) # Suse
	SUSE= true
	CC=/usr/lib64/mpi/gcc/openmpi/bin/mpicxx
else # Ubuntu etc.
	SUSE= false
	CC=/usr/bin/mpicxx
	LIB= -L/usr/lib/openmpi/lib
	INC= -I/usr/include/eigen3
endif

SDIR = src
ODIR = obj
PROG = bin/diagonalize
_OBJS = main.o MpiEnvironment.o MpiLoadbalancing.o MpiVariables.o MpiLoadbalancingSimple.o Numerov.o QuantumDefect.o DipoleMatrix.o SQLite.o\

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

all:	$(PROG)

$(PROG): $(OBJS)
	$(CC) -o $@ $(LDFLAGS) $(LIB) $(OBJS)
	@echo "Linking complete!"

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) $< -o $@ $(CFLAGS) $(INC)
	@echo "Compiled "$<" successfully!"

.PHONY: clean

clean:
	\rm -rf $(OBJS) $(PROG)

test:
	mpiexec -n 5 bin/diagonalize; # --hostfile /nodes/active_ompi

check: all test
