# Keeps a few frequent commands
#
cleantemp = rm -rf build; rm -f *.c

.PHONY : clean all build tarc extern libs

DIRS = include
TAR  = sourcecode.tar.gz

all: clean extern build

extern: 
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

commit: 
	git commit
	git push

libs: 
	python getlibs.py
