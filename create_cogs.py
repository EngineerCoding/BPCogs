from collections import Counter, OrderedDict
from multiprocessing.pool import ThreadPool

import psycopg2
import sys
import os


_BLAST_COMMAND = ("blastp -outfmt 6 -query {org1}.fa -db {org2}.fa " +
                  "-evalue 1e-10 -num_threads {n_threads} " +
                  "| awk '{{print $1,$2}}' | sed 's/ /;/g'" +
                  " > {org1}__{org2}")


def handle_blast_arguments():
    """ Handles the arguments given on the command line
    related to the BLAST settings. The different arguments are:
    --no-blast: marks the algorithm to not execute the BLAST commands
    --blast-threads <int>: the amount of threads which
    the BLAST command should be executed with.
    --no-blast-dbs: marks the algorithm to not execute the makeblastdb
    commands.

    Returns:
        1. Boolean: execute BLAST commands
        2. Boolean: execute makeblastdb commands
        3. Int: amount of BLAST threads
    """
    if "--no-blast" in sys.argv:
        return False, False, -1
    threads = 3
    if "--blast-threads" in sys.argv:
        index = sys.argv.index("--blast-threads")
        if (index + 1 < len(sys.argv) and
                sys.argv[index + 1].strip().isnumeric()):
            threads = max(1, int(sys.argv[index + 1].strip()))
    return True, "--no-blast-dbs" not in sys.argv, threads


def write_direct_hits(org1, org2):
    """ Retrieves the first BLAST hit from each protein and writes it
    to another file. The BLAST files have the format {org1}__{org2}
    which contain all the BLAST results for all the proteins from
    organism 1. For each protein the first result
    (and thus the best result) is retrieved and written to
    the {org1}_direct_{org2} file.

    Parameters:
        org1 (string): The name of the first organism
        org2 (string): The name of the second organism
    Returns:
        -
    """
    with open("{}_direct_{}".format(org1, org2), "w") as out:
        with open("{}__{}".format(org1, org2)) as infile:
            prev_prot = ""
            for result_line in infile:
                prot1, prot2 = result_line.split(";")
                if prot1 != prev_prot:
                    prev_prot = prot1
                    out.write(result_line)


def blast_organisms(organisms):
    """ Executes the BLAST and makeblastdb commands based on the
    settings given on the command line. Those settings are handled by
    the function 'handle_blast_arguments'. Not dependent on those
    settings is that the directional hits from the BLAST results are
    written using the 'write_direct_hits' function. This means that
    BLAST results must be available whether a BLAST is executed or not.

    Parameters:
        organisms (list of strings): The organisms which are used in
            this COG program.
    Returns:
        -
    """
    do_blast, make_databases, threads = handle_blast_arguments()
    if make_databases:
        for organism in organisms:
            os.system("makeblastdb -in {}.fa -dbtype prot".format(organism))
    for org1 in organisms:
        for org2 in organisms:
            if org1 != org2:
                if do_blast:
                    os.system(_BLAST_COMMAND.format(org1=org1, org2=org2,
                                                    n_threads=threads))
                write_direct_hits(org1, org2)


def insert_protein(cursor, protein_id, organism_id, protein_name, protein_seq,
                   protein2id):
    """ Inserts a protein into the database and also maps a protein
    name to its ID number. This ID is used later in the general
    algorithm.

    Parameters:
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
        protein_id (int): The ID which the to be inserted protein
            should have.
        organism_id (int): The ID of the organism which contains this
            protein.
        protein_name (string): The name of the protein.
        protein_seq (string): The amino acid sequence.
        protein2id (dictionary): The dictionary which has keys as names
            (practically only the accession) and values as numeric ID's.
    Returns:
        1. Int: A new protein ID for the next protein.
    """
    if not protein_name:
        return protein_id
    args = (protein_id, protein_name, protein_seq, organism_id)
    cursor.execute("""
        INSERT INTO protein (protein_id, name, sequence, organism)
        VALUES (%s, %s, %s, %s)""", args)
    protein2id[protein_name.split()[0]] = str(protein_id)
    return protein_id + 1


def write_protein_ids(protein2id, organism):
    """ Replaces the accession numbers in direct hit files by their
    corresponding ID in the protein table. This happens per organism,
    all file names which have '_direct_' and the organism name in it
    will be replaced. Those files, are later used to insert directly
    in the database.

    Parameters:
        protein2id (dictionary): The dictionary which has keys as names
            (practically only the accession) and values as numeric ID's.
        organism (string): The name of the organism.
    Returns:
        -
    """
    for file in os.listdir(os.getcwd()):
        if organism in file and "_direct_" in file:
            replace_side = int(file.split("_direct_")[1] == organism)
            with open(file) as infile:
                with open(file + ".tmp", "w") as outfile:
                    for line in infile:
                        prots = line.strip().split(";")
                        prots[replace_side] = protein2id[prots[replace_side]]
                        outfile.write(";".join(prots) + "\n")
            os.remove(file)
            os.rename(file + ".tmp", file)


def fill_protein_table(cursor, organism, org_id, prot_id):
    """ Fills the protein table one organism at a time. This will
    collect the full name (only the '>' is stripped off) and
    corresponding sequence. After the inserts are complete, the
    'write_protein_ids' function is called.

    Parameters:
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
        organism (string): The name of the organism.
        org_id (int): The database ID of the organism.
        prot_id (int): The next ID of a protein.
    Returns:
        1. Int: A new protein ID for the next protein.
    """
    protein2id = dict()
    last_prot = ""
    last_seq = ""

    with open(organism + ".fa") as infile:
        for line in infile:
            if line[0] == ">":
                prot_id = insert_protein(cursor, prot_id, org_id, last_prot,
                                         last_seq, protein2id)
                last_prot = line[1:].strip()
                last_seq = ""
            else:
                last_seq += line.strip()
    prot_id = insert_protein(cursor, prot_id, org_id, last_prot, last_seq,
                             protein2id)
    write_protein_ids(protein2id, organism)
    return prot_id


def fill_bidirectional_hits(connection, cursor):
    """ Fills the directionalhit table using all the *_direct_* files
    inserted into a temporary table. Then a query is executed to get
    all bidirectional hits and the temporary table will be removed.

    Parameters:
        connection (psycopg2 connection class): The connection with the
            database.
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
    Returns:
        -
    """
    os.system("cat *_direct_* > total")
    with open("total") as infile:
        cursor.copy_from(infile, "temporaryhit", sep=";",
                         columns=("protein_a", "protein_b"))
    connection.commit()
    os.remove("total")
    cursor.execute("""
        INSERT INTO directionalhit (protein_a, protein_b)
        SELECT protein_a, protein_b FROM (
          SELECT a.protein_a, a.protein_b
          FROM temporaryhit a
            INNER JOIN temporaryhit b
              ON a.protein_a = b.protein_b AND a.protein_b = b.protein_a
        ) AS combination_set
        WHERE protein_a < protein_b""")
    cursor.execute("DROP TABLE temporaryhit CASCADE")
    connection.commit()


def fill_database(connection, cursor, organisms):
    """ The main function which fills the database tables except for
    the cog and multiplesequencealignment tables. The tables are filled
    in the following order:
      1. organism table
      2. protein table
      3. directionalhit table

    Parameters:
        connection (psycopg2 connection class): The connection with the
            database.
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
        organisms (list of string): The organisms which are used in
            this COG program.
    Returns:
        -
    """
    bulk_insert = []
    for organism in organisms:
        bulk_insert.append(
            "INSERT INTO organism (name) VALUES ('{}')".format(organism))
    cursor.execute(";".join(bulk_insert))
    connection.commit()
    current_id = 1
    for organism in organisms:
        current_id = fill_protein_table(
            cursor, organism, organisms.index(organism) + 1, current_id)
    fill_bidirectional_hits(connection, cursor)


def remove_twogs(twogs, remove_containing):
    """ Removes twogs which contain a protein ID specified in the list
    of remove_containing. This method does NOT rewrite the twogs file.

    Parameters:
        twogs (list): A list of tuples which contain a bidirectional
            hit using protein ID's.
        remove_containing (list of ints): A list of protein ID's.
    Returns:
        -
    """
    remove_indices = []
    for index in range(len(twogs)):
        for remove_this in remove_containing:
            if remove_this in twogs[index]:
                remove_indices.append(index)
                break
    remove_indices.sort(reverse=True)
    for index in remove_indices:
        del twogs[index]


def write_read_twogs(twogs=None, mode="w", organism_id1=-1, organism_id2=-1):
    """ This method will always write the given twogs to a file called
    "twogs" with the specified mode. When twogs are not available, the
    organism_id1 and organism_id2 must be available to read the twogs
    from a file called "full_twogs". This reading involves only
    reading the twogs which are known from given organism ID's.

    Parameters:
        twogs (list): A list of tuples containing bidirectional hits.
        mode (string): The mode to write to the "twogs" file. Only "a"
            and "w" are permitted.
        organism_id1 (int): The organism ID of organism 1.
        organism_id2 (int): The organism ID of organism 2.
    Returns:
        -
    """
    if twogs is None:
        twogs = []
        organism_code_tuple = (organism_id1, organism_id2)
        with open("full_twogs") as infile:
            for line in infile:
                org_prot1, org_prot2, org_id1, org_id2 = tuple(map(
                    int, line.strip().split(";")))
                if ((org_id1, org_id2) == organism_code_tuple
                    or (org_id2, org_id1) == organism_code_tuple):
                        twogs.append((org_prot1, org_prot2))
    with open("twogs", mode) as outfile:
        for twog in twogs:
            outfile.write(";".join(map(str, twog)) + "\n")


def get_cog_proteins(cursor):
    """ Retrieves the cog number along with the protein ID's contained
    in it.

    Parameters:
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
    Returns:
        1. Dict: A dictionary with keys as cog ID and as
        values a list of protein ID's.
    """
    cursor.execute("""
        SELECT cog, protein_id
        FROM protein
        WHERE cog IS NOT NULL
        ORDER BY cog ASC""")
    cog_map = OrderedDict()
    for cog, protein_id in cursor.fetchall():
        if cog not in cog_map:
            cog_map[cog] = []
        cog_map[cog].append(protein_id)
    return cog_map


def update_cogs(cursor, twogs):
    """ Finds per cog a best matching protein and adds it to that cog.

    Parameters:
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
        twogs (list): A list of tuples containing bidirectional hits.
    Returns:
        1. List: A list of integers which represent protein ID's
        which should not be used in the 'new_cogs' function.
    """
    cog_map = get_cog_proteins(cursor)
    skip_proteins = []
    delete_twogs = []
    for cog in cog_map:
        cog_proteins = cog_map[cog]
        protein_counter = Counter()
        for twog in twogs:
            for cog_protein in cog_proteins:
                if cog_protein in twog:
                    index = int(twog[1] != cog_protein)
                    protein_counter.update(twog[index:index + 1])
                    delete_twogs.append(cog_protein)
        if len(protein_counter):
            common_protein = protein_counter.most_common(1)[0][0]
            cursor.execute("UPDATE protein SET cog = %s " +
                           "WHERE protein_id = %s", (cog, common_protein))
            remove_twogs(twogs, delete_twogs + [common_protein])
            delete_twogs.clear()
            skip_proteins.append(common_protein)
    return skip_proteins


def read_organism_proteins(organism_id, skip_proteins):
    """ Retrieves a set of protein ID's which belong to the given
    organism ID. This reduces the amount of loops which 'new_cogs'
    has to do. Proteins which are not in the bidirectionalhits will
    never be participating in a COG.

    Parameters:
        organism_id (int): The database ID of the organism.
        skip_proteins (list of ints): A list of integers which should
            are already in cogs according to the function 'update_cogs'.
    Returns:
        1. Set: A set of protein ID's which are used in 'new_cogs'.
    """
    proteins = set()
    with open("full_twogs") as infile:
        for line in infile:
            prot1, prot2, id1, id2 = tuple(map(int, line.strip().split(";")))
            if id1 == organism_id and prot1 not in skip_proteins:
                proteins.add(prot1)
            if id2 == organism_id and prot2 not in skip_proteins:
                proteins.add(prot2)
    return proteins


def new_cogs(cursor, twogs, organism_id, skip_proteins):
    """ Finds new cogs. A new cog is defined as 3 proteins which have
    bidirectional hits with each other:
        A - B
         \ /
          C

    Parameters:
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
        twogs (list): A list of tuples containing bidirectional hits.
        organism_id (int): The database ID of the organism.
        skip_proteins (list of ints): A list of integers which should
            are already in cogs according to the function 'update_cogs'
    Returns:
        -
    """
    cursor.execute("SELECT MAX(cog_id) FROM cog")
    cog_id = (cursor.fetchone()[0] or 0) + 1
    for protein_id in read_organism_proteins(organism_id, skip_proteins):
        found_twogs = []
        for twog in twogs:
            if protein_id in twog:
                found_twogs.append(twog[int(twog[1] != protein_id)])
        hidden_twogs = set()
        for found_twogs_index in range(len(found_twogs)):
            for combine_with in found_twogs[found_twogs_index + 1:]:
                twog_combo = (found_twogs[found_twogs_index], combine_with)
                if twog_combo in twogs or twog_combo[::-1] in twogs:
                    hidden_twogs.update([*twog_combo, protein_id])
        if len(hidden_twogs):
            cursor.execute("INSERT INTO cog (cog_id) VALUES (%s)", (cog_id,))
            cursor.execute("UPDATE protein SET cog = %s WHERE protein_id IN "
                           + str(tuple(hidden_twogs)), (cog_id,))
            remove_twogs(twogs, hidden_twogs)
            cog_id += 1


def find_cogs(connection, cursor, organisms):
    """ The main algorithm to find and update cogs. First it will dump
    all bidirectional hits from the database in a file called
    "full_twogs" along with the organism ID's of those proteins.
    Then it will read the twogs of two organisms and when the more or
    equal than 3 organism are contained in the twogs, the following
    sequence will happen:
      1. Read the twogs from the twogs file
      2. Call update_cogs
      3. Call new_cogs
      4. Write the twogs to remove the removed twogs completely.

    Parameters:
        connection (psycopg2 connection class): The connection with the
            database.
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
        organisms (list of string): The organisms which are used in
            this COG program.
    Returns:
        -
    """
    if os.path.exists("twogs"):
        os.remove("twogs")
    cursor.execute("""SELECT protein_a, protein_b, a.organism, b.organism
                      FROM directionalhit
                        INNER JOIN protein a ON a.protein_id = protein_a
                        INNER JOIN protein b ON b.protein_id = protein_b""")
    with open("full_twogs", "w") as outfile:
        for line_data in cursor:
            outfile.write(";".join(map(str, line_data)) + "\n")
    for organism1_id in range(len(organisms)):
        for organism2_id in range(len(organisms)):
            if organism1_id > organism2_id:
                write_read_twogs(mode="a", organism_id1=organism1_id + 1,
                                 organism_id2=organism2_id + 1)
        if organism1_id > 1:
            twogs = []
            with open("twogs") as infile:
                for line in infile:
                    twogs.append(tuple(map(int, line.strip().split(";"))))
            new_cogs(cursor, twogs, organism1_id + 1,
                     update_cogs(cursor, twogs))
            write_read_twogs(twogs=twogs, mode="w")
            connection.commit()


def get_msa(file):
    """ Parses an alignment file to fit the complete alignment in a
    string instead of a file.

    Parameters:
        file (string): The name of the file to parse
    Returns:
        1. String: the complete multiple sequence alignment
    """
    msa = OrderedDict()
    skip_line = True
    start_index = None
    with open(file) as infile:
        for line in infile:
            if skip_line:
                skip_line = False
                continue
            if len(line) > 1:
                lines = line.split()
                if len(lines) == 2:
                    if lines[0] not in msa:
                        msa[lines[0]] = ""
                    msa[lines[0]] += lines[1]
                    start_index = line.find(lines[1])
                else:
                    if "consensus" not in msa:
                        msa["consensus"] = ""
                    msa["consensus"] += line[start_index:-1]
    msa_text = ""
    for accession, sequence in msa.items():
        msa_text += '{}{}{}\n'.format(
            accession, " " * (start_index - len(accession)), sequence)
    return msa_text


def multiple_sequence_alignment(name_sequence_tuple, cog):
    """ A function which executes a multiple sequence alignment. First
    the sequence will be written to a file which then can be used by
    ClustalW. When ClustalW ran, the function 'get_msa' will parse the
    alignment file.

    Parameters:
        name_sequence_tuple (list): A list of tuples which contain the
            name of the sequence (0) and the sequence itself (1).
        cog (int): The cog ID.
    Returns:
        1. Int: the cog ID.
        2. String: the complete multiple sequence alignment.
    """
    files = "msa_file{c}.fa tree_file{c}.dnd output{c}.aln".format(
        c=cog).split()
    with open(files[0], "w") as outfile:
        for name, sequence in name_sequence_tuple:
            outfile.write(">" + name + "\n")
            for i in range(0, len(sequence), 80):
                outfile.write(sequence[i:i + 80] + "\n")
    os.system("clustalw -infile={} -newtree={} -outfile={} -align -quiet"
              .format(*files))
    msa_string = get_msa(files[2])
    for file in files:
        os.remove(file)
    return cog, msa_string


def do_multiple_sequence_alignments(cursor):
    """ The main function which can execute multiple sequence
    alignments by calling 'multiple_sequence_alignment'. This can be
    done in parallel (to increase the speed of the runtime) by adding
    '--msa-threads <int>' on the command line by replacing '<int>' by
    a positive integer greater than 0.
    The multiple sequence alignments from a cog are finally inserted in
    the database.

    Parameters:
        cursor (psycopg2 cursor class): The cursor to execute queries
            with.
    Returns:
        -
    """
    arguments = dict()
    cursor.execute("""
        SELECT name, sequence, cog
        FROM protein
        WHERE cog IS NOT NULL
        ORDER BY cog ASC""")
    for name, sequence, cog in cursor:
        if cog not in arguments:
            arguments[cog] = ([], cog)
        arguments[cog][0].append((name, sequence))
    msa_threads = 1
    if '--msa-threads' in sys.argv:
        index = sys.argv.index('--msa-threads')
        if (index + 1 < len(sys.argv) and
                sys.argv[index + 1].strip().isnumeric()):
            msa_threads = max(1, int(sys.argv[index + 1].strip()))
    pool = ThreadPool(msa_threads)
    msa_vars = pool.starmap(multiple_sequence_alignment, arguments.values())
    cursor.executemany("INSERT INTO multiplesequencealignment (cog, msa) " +
                       "VALUES (%s, %s)", msa_vars)


def main():
    if not os.path.exists("Project_files/Organismen.txt"):
        raise Exception("Project_files/Organismen.txt does not exist!")
    with open("Project_files/Organismen.txt") as f:
        organism_list = [organism.strip() for organism in f]
    os.chdir("Project_files")
    blast_organisms(organism_list)
    connection = psycopg2.connect(host="localhost", dbname="postgres",
                                  user="postgres", password="Password")
    cursor = connection.cursor()
    with open("../create_tables.sql") as sql:
        cursor.execute(sql.read())
    connection.commit()
    fill_database(connection, cursor, organism_list)
    find_cogs(connection, cursor, organism_list)
    do_multiple_sequence_alignments(cursor)
    connection.commit()
    cursor.close()
    connection.close()

if __name__ == "__main__":
    main()
