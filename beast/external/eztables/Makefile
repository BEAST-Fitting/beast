# Keeps a few frequent commands
#
PROJECT = 'tapy'
cleantemp = rm -rf build; rm -f *.c
TAR = $(PROJECT).tar.gz
DIRS= core backend

.PHONY : clean all build

all: clean build

build:  
	#for d in $(DIRS); do (cd $$d; $(MAKE) build );done
	python setup.py build
	$(cleantemp)


clean: 
	#for d in $(DIRS); do (cd $$d; $(MAKE) clean );done
	$(cleantmp)
	find . -name '*pyc' -exec rm -f {} \;

tar: clean
	
	tar zcf $(TAR) * --exclude=rabe
