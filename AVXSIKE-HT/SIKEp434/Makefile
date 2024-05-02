CC = gcc-9
FILES	=	$(wildcard src/*.c) $(wildcard src/*.S)
SOURCES = src/keccak_8x1w/KeccakP-1600-times8-SIMD512.o
OUTPUT = sike
FLAGS = -O2 -mavx512ifma -funroll-loops -fomit-frame-pointer -m64 -mbmi2 -march=native -fwrapv -mtune=native

all: src/keccak_8x1w/KeccakP-1600-times8-SIMD512.o
	$(CC) $(SOURCES) $(FILES) -o $(OUTPUT) $(FLAGS)

src/keccak_8x1w/KeccakP-1600-times8-SIMD512.o: src/keccak_8x1w/KeccakP-1600-times8-SIMD512.c \
	src/keccak_8x1w/align.h \
	src/keccak_8x1w/brg_endian.h \
	src/keccak_8x1w/KeccakP-1600-times8-SIMD512.c \
	src/keccak_8x1w/KeccakP-1600-times8-SnP.h \
	src/keccak_8x1w/KeccakP-1600-unrolling.macros \
	src/keccak_8x1w/SIMD512-config.h
	$(CC) $(FLAGS) -c src/keccak_8x1w/KeccakP-1600-times8-SIMD512.c -o $@

clean:
	rm src/keccak_8x1w/KeccakP-1600-times8-SIMD512.o
	rm -f sike


