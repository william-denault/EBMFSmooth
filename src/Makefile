# Compile all shared objects for TMB models
all: fit_theta.so opt_theta.so

# Compile fit_theta.cpp into fit_theta.so
fit_theta.so: fit_theta.cpp
	Rscript --vanilla -e "TMB::compile('fit_theta.cpp', '-O0 -g')"

# Compile opt_theta.cpp into opt_theta.so
opt_theta.so: opt_theta.cpp
	Rscript --vanilla -e "TMB::compile('opt_theta.cpp', '-O0 -g')"

# Clean up generated files
clean:
	rm -f *.so *.o
