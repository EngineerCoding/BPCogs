from collections import Counter

import psycopg2
import sys
import os


_BLAST_COMMAND = ("blastp -outfmt 6 -query {org1}.fa -db {org2}.fa " +
                  "-max_target_seqs 1 -num_threads {n_threads} " +
                  "| awk '{{print $1,$2}}' | sed 's/ /;/g'" +
                  " > {org1}__{org2}")


def handle_blast_arguments():
    if "--no-blast" in sys.argv:
        return False, False, -1
    threads = 3
    if "--blast-threads" in sys.argv:
        index = sys.argv.index("--blast-threads")
        if index + 1 < len(sys.argv):
            try:
                threads = max(1, int(sys.argv[index + 1]))
            except ValueError:
                raise Exception("Invalid argument to --blast-threads")
    return True, "--no-blast-dbs" not in sys.argv, threads


def write_direct_hits(org1, org2):
    with open("{}_direct_{}".format(org1, org2), "w") as out:
        with open("{}__{}".format(org1, org2)) as infile:
            prev_prot = ""
            for result_line in infile:
                prot1, prot2 = result_line.split(";")
                if prot1 != prev_prot:
                    prev_prot = prot1
                    out.write(result_line)


def blast_organisms(organisms):
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
    if not protein_name:
        return protein_id
    args = (protein_id, protein_name, protein_seq, organism_id)
    cursor.execute("INSERT INTO protein (protein_id, name, sequence," +
                   " organism) VALUES (%s, %s, %s, %s)", args)
    protein2id[protein_name.split()[0]] = str(protein_id)
    return protein_id + 1


def write_protein_ids(protein2id, organism):
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
    os.system("cat *_direct_* > total")
    with open("total") as infile:
        cursor.copy_from(infile, "temporaryhit", sep=";",
                         columns=("protein_a", "protein_b"))
    connection.commit()
    os.remove("total")
    cursor.execute(
        "INSERT INTO directionalhit (protein_a, protein_b) SELECT protein_" +
        "a, protein_b FROM (SELECT a.protein_a, a.protein_b FROM temporaryh" +
        "it a INNER JOIN temporaryhit b ON a.protein_a = b.protein_b AND a." +
        "protein_b = b.protein_a) AS combination_set WHERE protein_a < prot" +
        "ein_b")
    cursor.execute("DROP TABLE temporaryhit CASCADE")
    connection.commit()


def fill_database(connection, cursor, organisms):
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
    remove_indices = []
    for index in range(len(twogs)):
        for remove_this in remove_containing:
            if remove_this in twogs[index]:
                remove_indices.append(index)
                break
    remove_indices.sort(reverse=True)
    for index in remove_indices:
        del twogs[index]


def write_twogs(twogs=None, mode="w", organism_id1=-1, organism_id2=-1):
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


def update_cogs(cursor, twogs):
    cursor.execute("SELECT cog, protein_id FROM protein WHERE cog IS NOT NU" +
                   "LL ORDER BY cog ASC")
    cog_map = dict()
    for cog, protein_id in cursor.fetchall():
        if cog not in cog_map:
            cog_map[cog] = []
        cog_map[cog].append(protein_id)
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


def read_organism_proteins(organism_code, skip_proteins):
    proteins = set()
    with open("full_twogs") as infile:
        for line in infile:
            prot1, prot2, id1, id2 = tuple(map(int, line.strip().split(";")))
            if id1 == organism_code and prot1 not in skip_proteins:
                proteins.add(prot1)
            if id2 == organism_code and prot2 not in skip_proteins:
                proteins.add(prot2)
    return proteins


def new_cogs(cursor, twogs, organism_id, skip_proteins):
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
                write_twogs(mode="a", organism_id1=organism1_id + 1,
                            organism_id2=organism2_id + 1)
        if organism1_id > 1:
            twogs = []
            with open("twogs") as infile:
                for line in infile:
                    twogs.append(tuple(map(int, line.strip().split(";"))))
            new_cogs(cursor, twogs, organism1_id + 1,
                     update_cogs(cursor, twogs))
            write_twogs(twogs=twogs, mode="w")
            connection.commit()


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
    cursor.close()
    connection.close()

if __name__ == "__main__":
    main()
