#
# pyPEGASE MAKEFILE
# 
# Keeps a few frequent commands
#
.PHONY : clean all test tarc extern help
DIRS = include
TAR  = sourcecode.tar.gz

all: clean

build: 
	python setup.py build

install: 
	python setup.py install

clean : 
	find . -name '*pyc' -exec rm -f {} \;

.SILENT:
help:
	echo "\nOptions"
	echo "\t build   : python build package"
	echo "\t install : python install package"
	echo "\t clean   : clean repo from pyc so files"
	echo ""
		
