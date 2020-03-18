Strategy for two-step parallelisation:


1) Main parallelisation is for chains within a single model.

2) Secondary parallelisation ('poor man's parallelisation) is over species - this happens by a loop in the shell script sending variables
(species names) to multiple calls to qsub

3) instances of qsub take a species name as an environment variable and pass it to the call to R CMD BATCH

4) R script receives the variable using the commandArgs() option in R

5) At the moment job specific Rout files are created, but not unique job names



Note that if you create a new bash script it needs to be made executable:
# 

[olipes@cirrus scripts]$ chmod +x rr_FINAL.sh # make script executable

[olipes@cirrus scripts]$ ./rr_FINAL.sh # run script



Oli Pescott, 22/07/2016