all: ctmc_dive.so 

ctmc_dive.so: ctmc_dive.cpp
	Rscript --vanilla -e "TMB::compile('ctmc_dive.cpp','-O3 -g', safebounds=FALSE)"
	
clean:
	rm -rf *o 

