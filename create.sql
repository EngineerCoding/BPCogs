 CREATE TABLE organism
  (
     organism_id SERIAL NOT NULL,
     name        TEXT NOT NULL,
     CONSTRAINT organism_idpk PRIMARY KEY(organism_id)
  );

CREATE TABLE cog
  (
     cog_id SERIAL NOT NULL,
     msa    TEXT NOT NULL,
     CONSTRAINT cog_idpk PRIMARY KEY (cog_id)
  );

CREATE TABLE protein
  (
     protein_id SERIAL NOT NULL,
     name       TEXT NOT NULL,
     sequence   TEXT NOT NULL,
     organism   SERIAL NOT NULL,
     cog        INTEGER,
     CONSTRAINT protein_idpk PRIMARY KEY(protein_id),
     CONSTRAINT organism_fk FOREIGN KEY(organism) REFERENCES organism(organism_id),
     CONSTRAINT cog_fk FOREIGN KEY(cog) REFERENCES cog(cog_id)
  );

CREATE TABLE directionalhit
  (
     hit_id    SERIAL NOT NULL,
     protein_a INTEGER NOT NULL,
     protein_b INTEGER NOT NULL,
     CONSTRAINT hit_idpk PRIMARY KEY(hit_id),
     CONSTRAINT protein_afk FOREIGN KEY(protein_a) REFERENCES protein(protein_id),
     CONSTRAINT protein_bfk FOREIGN KEY(protein_b) REFERENCES protein(protein_id)
  );  
