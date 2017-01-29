import psycopg2
import os
conn = psycopg2.connect(host='localhost',dbname='postgres',user='postgres',password='hoi')
cursor = conn.cursor()



sql_organism = 'CREATE TABLE Organism (organism_ID SERIAL \
NOT NULL, Naam TEXT NOT NULL, CONSTRAINT Organism_IDPK\
 PRIMARY KEY(Organism_ID))'
cursor.execute(sql_organism)

sql_cog = 'CREATE TABLE cog (cog_ID SERIAL NOT NULL, msa TEXT\
 NOT NULL, CONSTRAINT cog_IDPK PRIMARY KEY (cog_ID))'
cursor.execute(sql_cog)

sql_protein = 'CREATE TABLE protein (protein_ID SERIAL NOT NULL, name TEXT NOT NULL, sequence TEXT NOT NULL, organism SERIAL NOT NULL, cog SERIAL NOT NULL, CONSTRAINT protein_IDPK PRIMARY KEY(protein_ID),CONSTRAINT organism_FK FOREIGN KEY(organism) REFERENCES organism(organism_ID),CONSTRAINT cog_FK FOREIGN KEY(cog) REFERENCES cog(cog_ID))'
cursor.execute(sql_protein)

sql_directionalhit = 'CREATE TABLE directionalhit (hit_ID SERIAL NOT NULL, protein_A INTEGER NOT NULL, protein_B INTEGER NOT NULL, CONSTRAINT hit_IDPK PRIMARY KEY(hit_ID), CONSTRAINT protein_AFK FOREIGN KEY(protein_A) REFERENCES protein(protein_ID), CONSTRAINT protein_BFK FOREIGN KEY(protein_b) REFERENCES protein(protein_ID))'
cursor.execute(sql_directionalhit)

conn.commit()

####" 

####
