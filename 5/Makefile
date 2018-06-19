CC=g++
OPT= -std=c++11 -O3 -I/usr/include/eigen3 -fopenmp
blatt_nummer=5
zip_inhalt=2.cpp 2.py build/main.pdf

all: build/2
latex: all build/main.pdf

build/2: 2.cpp | build
	${CC} ${OPT} 2.cpp -o build/2
	build/2
	python 2.py

include ../Makefile_template.mk
