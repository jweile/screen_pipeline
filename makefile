R=/home/jweile/bin/Rscript
NOW := $(shell date +"%y-%m-%d")
STABLE=v0.1
PROJECT=screen_pipeline

#default target: unit test and build
all: build_stable

# #run all unit tests
# test: setup_test
# 	cd testbench;\
# 	for test in *_test.R; do \
# 		$(R) $$test >$${test}.log;\
# 	done

# #setup a test directory
# setup_test:
# 	mkdir -p testbench/lib/
# 	mkdir -p testbench/res/
# 	mkdir -p testbench/html/
# 	cp src/main/R/* testbench/lib/
# 	cp src/main/resources/* testbench/res/
# 	cp src/main/html/* testbench/html
# 	cp src/test/R/* testbench/
# 	cp src/test/bash/* testbench/
# 	cp src/test/resources/* testbench/res/

# #delete test data
# clean:
# 	rm -r testbench/

#load the latest version of the code
latest:
	hg update

#load the last stable version of the code
stable:
	hg update $(STABLE)

#build a zip file with the final product
_build: 
	#create build directory
	mkdir -p $(PROJECT)/lib/
	# mkdir -p $(PROJECT)/html/
	mkdir -p $(PROJECT)/res/
	#copy code into dir structure
	cp src/main/R/* $(PROJECT)/lib/
	cp src/main/resources/* $(PROJECT)/res/
	# cp src/main/html/* $(PROJECT)/html/
	cp src/main/bash/* $(PROJECT)/
	#interpolate binary locations
	src/make/interpolate.sh bins.cfg $(PROJECT)/
	#setup permissions
	chmod u+x $(PROJECT)/lib/*.R
	chmod u+x $(PROJECT)/*.sh
	#zip build
	zip $(PROJECT)_$(NOW).zip -r $(PROJECT)/ 
	#delete build directory
	rm -r $(PROJECT)/

build_latest: latest _build

build_stable: stable _build

#install the built contents as well as required dependencies
install: 
	$(eval R=`grep Rscript bins.cfg|cut -d, -f2`)
	$(R) src/make/install_dependencies.R
	unzip $(PROJECT)_$(NOW).zip -d $${HOME}
	# mkdir -p $${HOME}/www/html/
	# if ! [ -h $${HOME}/www/html/$(PROJECT) ]; then \
	# 	ln -s $${HOME}/$(PROJECT)/ $${HOME}/www/html/;\
	# fi
