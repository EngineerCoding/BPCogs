
with open("Bpcogs/all_proteins.txt", 'r') as all_proteins:
    protein_list = [[] for i in range(0, 10)]
    org_count = 0

    for protein in all_proteins:
        if protein[0] == '#':
            protein_list[org_count].append(protein.strip('# \n'))
            org_count += 1
            continue
        if protein[0] == '>':
            protein_list[org_count - 1].append(protein.strip())
