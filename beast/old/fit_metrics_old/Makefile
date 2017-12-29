# Keeps a few frequent commands
#
cleantemp = rm -rf build; #rm -f *.c

.PHONY : clean all build

all: clean build

build:  
	python setup.py build_ext --inplace
	$(cleantemp)

build3:
	python3 setup.py build_ext --inplace
	$(cleantemp)

clean: 
	$(cleantmp)
	find . -name '*pyc' -exec rm -f {} \;
	rm -f *.so
	rm -rf __pycache__

