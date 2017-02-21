import psycopg2
conn = psycopg2.connect(host='localhost',dbname='postgres',user='postgres',password='hoi')
cursor = conn.cursor()

sql_executor = open('create.sql','r')
cursor.execute(sql_executor.read())
sql_executor.close()

with open('Bpcogs/all_proteins.txt', 'r') as all_proteins:
    protein_list = [[] for i in range(0, 10)]
    org_count = 0

    for protein in all_proteins:
        if protein[0] == '#':
            protein_list[org_count].append(protein.strip('# \n'))
            org_count += 1
            continue
        if protein[0] == '>':
            protein_list[org_count - 1].append(protein.strip())


for organism in protein_list:
    query = 'INSERT INTO organism(name) VALUES (organism[organism][0])'
    cursor.execute(query)

for x in range(0, len(protein_list)):
    for y in range(0, len(protein_list)):
        insert_protein = 'INSERT INTO protein(name, organism) VALUES (protein_list[x][y], protein_list[x][0]'
        cursor.execute(insert_protein)

with open('Bpcogs/bidirectinalhits.txt', 'r') as bdh:
    for line in bdh:
        if line[0] == '>':
            continue
        protein = line.rstrip().split(';')
        insert_bdh = 'INSERT INTO directionalhit(protein_a, protein_b) VALUES (protein[0], protein[1])'
        cursor.execute(insert_bdh)

conn.commit()

####" 

####
