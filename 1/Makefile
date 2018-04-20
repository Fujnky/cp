CC=g++
OPT= -std=c++11 -O3
blatt_nummer=1
zip_inhalt=1.cpp 1.py 2.cpp build/main.pdf Makefile

all: build/1 build/2
latex: all build/main.pdf

build/1: 1.cpp 1.py | build
	@printf "Aufgabe 1:\n"
	${CC} ${OPT} 1.cpp -o build/1
	@build/1
	python3 1.py

build/2: 2.cpp | build
	@printf "\n\nAufgabe 2:\n"
	${CC} ${OPT} 2.cpp -o build/2
	@build/2

zip: all latex | build
	mkdir -p build/zipfile
	cp -R ${zip_inhalt} build/zipfile
	cd build/zipfile && zip -r -9 ../Blatt$(blatt_nummer)_Dag-Björn\ Hering_Lars\ Funke.zip .
	rm -r build/zipfile

build/main.pdf: FORCE | build
	  TEXINPUTS="$(call translate,build:)" \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --xelatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
