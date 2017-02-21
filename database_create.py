import psycopg2
import os
conn = psycopg2.connect(host='localhost',dbname='postgres',user='postgres',password='hoi')
cursor = conn.cursor()

sql_executor = open ('create.sql','r')
cursor.execute(sql_executor.read())
sql_executor.close()

organismen_txt = open ('Organismen.txt','r')
for organisme in organismen_txt:
    'INSERT INTO Enzym_05(Naam) VALUES (organisme)'


conn.commit()

####" 

####
