# Keeps a few frequent commands
#
cleantemp = rm -rf build; rm -f *.c

.PHONY : clean all build tarc extern libs

DIRS = include
TAR  = sourcecode.tar.gz


extern: build
	for d in $(DIRS); do (cd $$d; $(MAKE) build );done

clean: 
	find . -name '*pyc' -exec rm -f {} \;
	find . -name '*so' -exec rm -f {} \;
	rm -f $(TAR)
	for d in $(DIRS); do (cd $$d; $(MAKE) clean );done

tar: clean $(TAR)

$(TAR):
	tar zcf $(TAR) * --exclude=libs

build:  
	python setup.py build_ext --inplace
	$(cleantemp)

libs: 
	python getlibs.py

all: clean extern build libs
