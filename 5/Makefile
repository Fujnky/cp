CC=g++
OPT= -std=c++11 -O3 -larmadillo -I/usr/include/eigen3 -fopenmp
blatt_nummer=4
zip_inhalt=1.cpp 1.py 1_beta.py build/main.pdf build/anim_q2.mp4 build/anim_q3.mp4

all: build/1
latex: all build/main.pdf

build/1: 1.cpp | build
	${CC} ${OPT} 1.cpp -o build/1


#build/2: 2.cpp | build
#	@printf "Aufgabe 2:\n"
#	${CC} ${OPT} 2.cpp -o build/2
#	build/2
#	python3 2.py

include ../Makefile_template.mk
