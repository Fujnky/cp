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

zip: all | build
	mkdir -p build/zipfile
	cp -R $(zip_inhalt) build/zipfile
	sed '$$d' Makefile > build/zipfile/Makefile
	cat ../Makefile_template.mk >> build/zipfile/Makefile
	cd build/zipfile && zip -r -9 ../Blatt$(blatt_nummer)_Dag-Björn\ Hering_Lars\ Funke.zip .

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean zip
