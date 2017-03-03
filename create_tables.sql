DROP TABLE IF EXISTS organism CASCADE;
CREATE TABLE organism
  (
     organism_id SERIAL NOT NULL,
     name        TEXT NOT NULL,
     CONSTRAINT organism_idpk PRIMARY KEY(organism_id)
  );

DROP TABLE IF EXISTS cog CASCADE;
CREATE TABLE cog
  (
     cog_id SERIAL NOT NULL,
     CONSTRAINT cog_idpk PRIMARY KEY (cog_id)
  );

DROP TABLE IF EXISTS protein CASCADE;
CREATE TABLE protein
  (
     protein_id SERIAL NOT NULL,
     name       TEXT NOT NULL,
     sequence   TEXT NOT NULL,
     organism   INTEGER NOT NULL,
     cog        INTEGER,
     CONSTRAINT protein_idpk PRIMARY KEY(protein_id),
     CONSTRAINT organism_fk FOREIGN KEY(organism) REFERENCES organism(organism_id),
     CONSTRAINT cog_fk FOREIGN KEY(cog) REFERENCES cog(cog_id)
  );

DROP TABLE IF EXISTS temporaryhit CASCADE;
CREATE TABLE temporaryhit
  (
      protein_a INTEGER NOT NULL,
      protein_b INTEGER NOT NULL
  );

DROP TABLE IF EXISTS directionalhit CASCADE;
CREATE TABLE directionalhit
  (
     hit_id    SERIAL NOT NULL,
     protein_a INTEGER NOT NULL,
     protein_b INTEGER NOT NULL,
     CONSTRAINT hit_idpk PRIMARY KEY(hit_id),
     CONSTRAINT protein_afk FOREIGN KEY(protein_a) REFERENCES protein(protein_id),
     CONSTRAINT protein_bfk FOREIGN KEY(protein_b) REFERENCES protein(protein_id)
  );

DROP TABLE IF EXISTS multiplesequencealignment CASCADE;
 CREATE TABLE multiplesequencealignment
  (
     msa_id    SERIAL NOT NULL,
     msa       TEXT NOT NULL,
     cog       INTEGER,
     CONSTRAINT msa_idpk PRIMARY KEY(msa_id),
     CONSTRAINT cog_fk FOREIGN KEY(cog) REFERENCES cog(cog_id)
  );
