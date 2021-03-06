# BOVERI-466: Generate Strelka indel calls for all MiSeq v4.0 runs
Needed for BOVERI-462: Generate metadata necessary for indel testing

A list of MiSeq runs for amplicons manifest v4.0 was obtained from <a href="https://docs.google.com/spreadsheets/d/1vIIKf5DTvHQy_mG7rCq0oCy-z82kj7bmv88wYcsWWCo/edit#gid=779463918">RUN_SUMMARIES</a>, on December 17, 2020, up until CG001Qv42Run248, and stored in ```data/RUN_SUMMARIES_MiSeq_v40.csv```. A list of excluded runs and samples was obtained from

using <a href="https://github.com/contextual-genomics/Bioinformatics/blob/dev/longitudinal_analysis/cohort_details/master_RD_excluded_RUNS.csv">master_RD_excluded_RUNS.csv</a> and <a href="https://github.com/contextual-genomics/Bioinformatics/blob/dev/longitudinal_analysis/cohort_details/master_RD_excluded_samples.csv">master_RD_excluded_samples.csv</a>s. This resulted in the files
```data/EXCLUDED_RUNS_v40.csv``` and ```data/EXCLUDED_SAMPLES_v40.csv```.


For each non-excluded MiSeq run, the latest_run information was obtained using the script <a href="https://github.com/contextual-genomics/Bioinformatics/blob/master/Operations/get_pipeline_latest_output.py">get_pipeline_latest_output.py</a>. This was done through the command  
```python bin/main.py parameters.tsv log/latest_runs_v40.log latest_runs > log/latest_runs_v40.out 2>&1```  
The result is in the file ```out/RUN_SUMMARIES_MiSeq_v40_latest_run.csv```.

The only run for which I could not get the results because the dump_variants did not work is CG001Qv42Run223. 

The variants dump files were obtained using a fixed version of the script <a href="https://github.com/contextual-genomics/biosys/blob/rd/rd_analysis/dump_variants.py">dump_variants.py</a>.  This was done through the command     
```python bin/main.py parameters.tsv log/dump_files_v40.log dump_files > log/dump_files_v40.out 2>&1```  
The TSV file containing the list of paths dump files indexed by run is ```out/dump_files.tsv```.  
Due to the fact the current operational dump_variants.py sxript contains a bug related to the coverage computation, I ran a local, fixed, version and stored the dump files locally, in ```out/dumpvariants_output```, that was not commited to github due to its size.



Finally the main variants CSV file, that contains information about all the indels observed in at least one of the non-excluded samples from a run in a matching run pair was generated by the command  
```python bin/main.py parameters.tsv log/variants_v40.log variants > log/variants_v40.out 2>&1```  
It resulted in a CSV file, ```out/variants_v40.csv```
