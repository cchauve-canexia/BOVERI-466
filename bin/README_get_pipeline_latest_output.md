# How to use get_pipeline_latest_output.py
The script provides you with the latest pipeline version which is tagged as IS_LATEST_RUN given the run, like 191104_M03829_0264_000000000-CN82P. 

## Input:
CSV file with "Raw Data" column that list the runs.

## Output 
CSV file with "Raw Data" and "latest_run" column that list the runs and their correponding *latest_pipeline_run* 

If there is there is no piepline run or there is more than one tagged with "IS_LATEST_RUN", you will get WARNINGS.
