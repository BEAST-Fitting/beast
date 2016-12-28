# Keeps a few frequent commands
#
cleantemp = rm -rf build; rm -f interp.c

.PHONY : clean all build

all: clean build

build:  
	python setup.py build_ext --inplace
	$(cleantemp)


clean: 
	$(cleantmp)
	find . -name '*pyc' -exec rm -f {} \;
	rm -f interp.c *.so

