# BPCogs

This file is meant to be viewed as a Markdown file, as Github displays it. If you are currently not reading
this file on Github, I would recommend to do so: https://github.com/EngineerCoding/BPCogs-2016-2017/blob/master/README.md.
Also, in the same directory of this file a PDF version can be found, so no internet has to be wasted to watch this file.


The script is based on the 'cogs_builder.py' file as programmed by an intern, but our script, 'create_cogs.py',
will do everything from scratch. The cogs_builder.py script will only find the cogs assuming all the files with
best directional hits are available.


The create_cogs.py script will basically do 4 things:
1. BLAST all proteomes against each other.
2. Find the best directional hits from the BLAST results.
3. Load the following data into a PostgreSQL database:
    1. All the organisms with a name along with an identifier
    2. All proteins, which contain the name, sequence and organism from which it originates from.
    3. A table with bidirectional hits
    4. A table which contain the cog identifiers. The identifiers are set on the proteins which belong
        in that cog. (not done in this step)
    5. And last but not least, a table with all multiple sequence alignments from the cogs. (not done in this step)
4. Find the cogs
5. Do the multiple sequence alignments


Since BLASTing the proteomes can take quite a while, it can be skipped and other additional arguments can be set
on the command line:

| Command line argument | Additional parameter | Description |
|-|-|-|
|--no-blast| - | Do not BLAST the proteomes, the results should already be available (no extra check is done for this!)|
|--blast-threads | int: amount of threads | The amount of threads BLAST should be using |
|--no-blast-dbs | - | Does not create database from the fasta files, those should then already available (no extra check is done for this)|
|--msa-threads | int: amount of threads | The amount of threads which should be used to create multiple sequence alignments with |


The required files to run the script with the current 'Organismen.txt', can be found in the in file 'Proteomes.zip'. This file should be 
extracted to the Project_files directory. When one also wants to skip the BLAST, the contents of 'blast_results.zip' also should be
copied to the Project_files directory.


The recommended way to run the program, given that BLAST results are already present, is using the following command:
```
python3 create_cogs.py --no-blast --msa-threads 4
```


It is possible that the script will result in a error, which is caused by the psycopg2 module. Most likely the credentials of the
database are incorrect, which need to be edited to have the correct credentials. This is done in the main function of the script.
Currently the option does not exist to give credentials to the program using the command line.


The finding cogs part of this script, follows the same algorithm as the cogs_builder.py with a few improvements:
* CogUpdate (update_cogs in the new script):
    1. The cogs are fetched along with the proteins which are in that cog, reducing search time. No actual test
    is done to compare those times, but theoretically the current implementation is faster due to the use of a
    database.
* NewCogFind (new_cogs in the new script):
    1. The list of proteins which gets iterated through, is optimized since only the proteins which are in the
    bidirectional hits are used instead of the full proteomes of all organisms.
    2. The HiddenTwogs loop is also optimised by checking in each iteration for (a, b) and (b, a) if those
    are in twogs, reducing the amount of looping which needs to be done. The loop is almost cut in half
    compared to the function in the original script.
