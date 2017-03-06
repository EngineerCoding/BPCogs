import psycopg2
import sys
import os


_BLAST_COMMAND = ("blastp -outfmt 6 -query {org1} -db {org2} " +
                  "-evalue 0.0000000001 -num_threads {n_threads} " +
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
    if not do_blast:
        return
    if make_databases:
        for organism in organisms:
            os.system("makeblastdb -in {}.fa -dbtype prot".format(organism))
    for org1 in organisms:
        for org2 in organisms:
            if org1 != org2:
                os.system(_BLAST_COMMAND.format(org1=org1, org2=org2,
                                                n_threads=threads))
                write_direct_hits(org1, org2)


def insert_protein(cursor, protein_id, organism_id, protein_name, protein_seq):
    if not protein_name:
        return protein_id
    args = (protein_id, organism_id, protein_name, protein_seq)
    cursor.execute("INSERT INTO protein (protein_id, name, sequence," +
                   " organism) VALUES (%s, %s, %s, %s)", args)
    return protein_id + 1


def write_protein_ids(protein2id, organism):
    for relative in os.listdir(os.getcwd()):
        if (os.path.isfile(relative) and organism in relative and
                    "_direct_" in relative):
            replace_side = int(relative.split("_direct_")[1] == organism)
            with open(relative + ".ids", "w") as out:
                with open(relative) as infile:
                    for line in infile:
                        prots = line.strip().split(";")
                        prots[replace_side] = protein2id[prots[replace_side]]
                        out.write(";".join(prots) + "\n")


def fill_protein_table(cursor, organism, org_id, prot_id):
    protein2id = dict()
    last_prot = ""
    last_seq = ""

    with open(organism + ".fa") as infile:
        for line in infile:
            if line[0] == ">":
                prot_id = insert_protein(cursor, prot_id, org_id, last_prot,
                                         last_seq)
                last_prot = line[1:].strip()
                last_seq = ""
            else:
                last_seq += line.strip()
    prot_id = insert_protein(cursor, prot_id, org_id, last_prot, last_seq)
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
    connection.commit()
    cursor.execute("DROP TABLE temporaryhit CASCADE")


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


def find_cogs(connection, cursor):
    pass


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
    find_cogs(connection, cursor)
    cursor.close()
    connection.close()

if __name__ == "__main__":
    main()
