CC=g++
OPT= -std=c++11 -O3
blatt_nummer=3a
zip_inhalt=1.cpp 1.py build/main.pdf

all: build/1# build/2
latex: all build/main.pdf

build/1: 1.cpp 1.py | build
	@printf "Aufgabe 1:\n"
	${CC} ${OPT} 1.cpp -o build/1
	@build/1
	python3 1.py

#build/2: 2.cpp 2.py | build
#	@printf "\n\nAufgabe 2:\n"
#	${CC} ${OPT} 2.cpp -o build/2
#	@build/2
#	python3 2.py

include ../Makefile_template.mk
