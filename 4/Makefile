CC=g++
OPT= -std=c++11 -O3 -larmadillo -fopenmp
blatt_nummer=4
zip_inhalt=1.cpp 1.py 1_beta.py build/main.pdf build/anim_q2.mp4 build/anim_q3.mp4

all: build/1 build/beta_q2.pdf build/q3_t0.pdf
latex: all build/main.pdf

build/1: 1.cpp | build
	${CC} ${OPT} 1.cpp -o build/1

build/beta_q2.pdf: build/1 1_beta.py | build
	build/1 --beta 2 | python 1_beta.py
	build/1 --beta 3 | python 1_beta.py
	build/1 --beta 4 | python 1_beta.py

build/q3_t0.pdf: build/1 1.py | build
	build/1 --history 2 | python 1.py
	build/1 --history 3 | python 1.py

#build/2: 2.cpp | build
#	@printf "Aufgabe 2:\n"
#	${CC} ${OPT} 2.cpp -o build/2
#	build/2
#	python3 2.py

include ../Makefile_template.mk
