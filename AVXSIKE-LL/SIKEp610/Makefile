CC = gcc-9
FILES	=	$(wildcard src/*.c) $(wildcard src/*.S)
OUTPUT = sike
FLAGS = -O2 -mavx512ifma -fomit-frame-pointer -m64 -mbmi2 -march=native -fwrapv -mtune=native

KAT_FILES = ./KAT/PQCtestKAT_kem610.c \
            ./KAT/rng/rng.c \
            ./KAT/aes/aes.c \
            ./KAT/aes/aes_c.c \
            ./src/utils.c \
            ./src/fips202.c \
            ./src/fp.c \
            ./src/fpx.c \
            ./src/curve.c \
            ./src/sidh.c \
            ./src/sike.c \
            ./src/fp_x64_asm.S
KAT_OUTPUT = kat

all: 
	$(CC) $(FILES) -o $(OUTPUT) $(FLAGS)

kat:
	$(CC) $(KAT_FILES) -o $(KAT_OUTPUT) $(FLAGS)

clean:
	rm -f sike kat
