CC=g++
OPT= -std=c++11 -O3 -I/usr/include/eigen3
blatt_nummer=7
zip_inhalt=1.cpp 1.py 2.py d.py build/main.pdf

all: build/1
latex: all build/main.pdf

build/1: 1.cpp 1.py | build
	${CC} ${OPT} 1.cpp -o build/1
	build/1 | python 1.py
include ../Makefile_template.mk
