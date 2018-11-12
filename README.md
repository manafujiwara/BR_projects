# Summary
Matlab codes used to analyse in a paper (Fujiwara et al., 2018, PLoS ONE),
which investigated use of Optokynetic nystagmus (a kind of eye movemnet) to
predict perception during binocular rivalry in Parkinson's disease patients.

## Discription
- The codes suppose to analyse data collected in an experiment using codes in
 https://github.com/manafujiwara/BR_projects/experiments

- Put all the analysis codes into a folder (leave nested codes, for example,
aux_files, as they are) that is going to be the home directory. For all
directories made automatically, please look at "setDir.m".

- First, Run "[DIR, subID] = startAnalysis.m" to obtain variables "DIR" and "subID",
which represent directories and subject ID. These variables are required to run
farther processes.

- Open "batch_BR_analysis.m" to start analysis, and choose which process to run.
