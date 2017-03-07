# FastQlyzer
!!!!!!!!!!!!SEE FASTQLYZER REPORT FILES!!!!!!!!!!!!!!!

Python tool that provides a way to do some quality and content control
checks on raw sequence data coming from high throughput sequencing pipelines.
It provides a modular set of analyses which you can use to give a quick 
impression of whether your data has any problems which you should be aware
before doing any further analysis.
FastQlyzer is wrapped on SBG platform and tested on public data.
Main functions of FastQlyzer are:
                  -Analyzing data by creating multiple graphs and tables that represent results of analyzing
                  -Export of results to a .pdf report file


Usage: FastQlyzer.py [-h] [--qt=quality_threshold] FASTQ1 FASTQ2

Options:
    -h --help
    --qt Quality threshold between bad and good quality. [Default: 9]
