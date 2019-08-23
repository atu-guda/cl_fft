CFLAGS=-O2 -g3 -Wall
CXXFLAGS=-O2 -g3 -Wall -std=c++17 -fgnu-keywords -Weffc++ -Woverloaded-virtual -lm -lfftw3

PGM=cl_fft
all: $(PGM)

clean:
	rm -f $(PGM) *.o core* *.dvi *.ps *.log

.PHONY: clean

