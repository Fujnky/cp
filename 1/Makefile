CC=g++
OPT= -std=c++11 -O3
blatt_nummer=1
zip_inhalt=1.cpp 1.py 2.cpp build/main.pdf

all: build/1 build/2

build/1: 1.cpp 1.py | build
	@printf "Aufgabe 1:\n"
	${CC} ${OPT} 1.cpp -o build/1
	build/1
	python3 1.py

build/2: 2.cpp | build
	@printf "\n\nAufgabe 2:\n"
	${CC} ${OPT} 2.cpp -o build/2
	build/2

include ../Makefile_template.mk
