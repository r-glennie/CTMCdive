all: ctmc_dive.so 

ctmc_dive.so: ctmc_dive.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('ctmc_dive.cpp','-O0 -g -Wno-ignored-attributes', safebounds=FALSE)"
	
clean:
	rm -rf *o 

