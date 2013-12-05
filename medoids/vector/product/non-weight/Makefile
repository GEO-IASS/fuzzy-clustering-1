CC=g++
CFLAGS=-c -Wall -O2 -g -W
LDFLAGS=
SOURCES=src/Main.cpp src/fmm/Algoritmo.cpp src/fmm/Array.cpp src/fmm/Cluster.cpp src/fmm/Dados.cpp src/fmm/Imprime.cpp src/fmm/Includes.cpp src/fmm/Repositorio.cpp src/fmm/Resultado.cpp src/fmm/Tabela.cpp src/sodas/MonteCarlo.cpp src/sodas/Pattern.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Run.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf src/*.o
	rm -rf src/fmm/*.o
	rm -rf src/sodas/*.o


