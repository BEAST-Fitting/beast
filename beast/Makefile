# Keeps a few frequent commands
#
cleantemp = rm -rf build; rm -f *.c

.PHONY : clean all build tarc extern libs help
.SILENT: help

DIRS = include proba
TAR  = sourcecode.tar.gz


help:
	echo "Make has many targets."
	echo "'all' will make all the code updates from the git repository, compile necessary files, and download libraries"
	echo "'all-but-libs' is equivalent to 'all' but does not download the libraries"

all-but-libs: clean gitmain gitdeps extern build

all: clean gitmain gitdeps extern build libs

gitdeps:
	git submodule init
	git submodule update --recursive
	git submodule foreach git pull origin master

gitmain:
	git pull

gitbuild: gitdeps gitmain


extern: build
	for d in $(DIRS); do (cd $$d; $(MAKE) build );done

clean: 
	$(cleantmp)
	find . -name '*pyc' -exec rm -f {} \;
	rm -f *.so
	rm -rf __pycache__
	for d in $(DIRS); do (cd $$d && $(MAKE) clean ); done
	rm -f $(TAR)

tar: clean $(TAR)

$(TAR):
	tar zcf $(TAR) * --exclude=libs

proba:  
	+$(MAKE) -C proba build

build: proba 
	python setup.py build_ext --inplace
	$(cleantemp)

libs: 
	python getlibs.py

