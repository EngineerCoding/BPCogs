# BPCogs

This file is meant to be viewed as Markdown file, as github displays it. If you are currently not reading
this file on github, I would recommend to do so: https://github.com/EngineerCoding/BPCogs-2016-2017/blob/master/README.md.


The script is based on cogs_builder.py file as built by an intern, but our script will do everyhting from
scratch. The cogs_builder.py script will only find the cogs assuming all the files with best directional
hits are available.


The create_cogs.py script will basically do 4 things:
1. BLAST all proteomes against each other.
2. Find the best directional hits from the BLAST results.
3. Load data into a PostgreSQL database:
    1. All the organisms with a name along with an identifier
    2. All proteins, which contain the name, sequence and organism from which it originates from.
    3. A table with bidirectional hits
    4. A table which contain the cog identifiers. The identifiers are set on the proteins which belong
        in that cog. (not done in this step)
    5. And last but not least, a table with alle multiple sequence alignments from the cogs. (not done in this step)
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



The finding cogs part of this script, follows the same algorithm as the cogs_builder.py with a few improvements:
* CogUpdate (update_cogs in the new script):
    1. The cogs are fetched along with the proteins which are in that cog, reducing search time. No actual test
    is done to compare those times, but theoretically the current implementation is faster due to the use of a
    database.
* NewCogFind (new_cogs in the new script):
    1. The list of proteins which gets iterated through, is optimised since only the proteins which are in the
    bidirectional hits are used instead of the full proteomes of all organisms.
    2. The VebrogenTwogs loop is also optimised by checking in each iteration for (a, b) and (b, a) if those
    are in twogs, reducing the amount of looping which needs to be done. The loop is almost cut in half
    compared by the original function.
