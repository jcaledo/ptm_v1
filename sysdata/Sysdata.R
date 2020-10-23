<<<<<<< HEAD
path1 <- "/Users/JCA/Dropbox/Investigacion/R_ptm/ptm_bis/sysdata"
path2 <- "/Users/JCA/Dropbox/Investigacion/R_ptm/ptm_bis/ptm/Rptm"

setwd(path1)

load(file = "./aai.Rda") # Amino acids properties
load(file = "./ac_db.Rda") # Acetylation DB 2019
load(file = "./d.Rda") # PDB-UniProt 2019
load(file = "./dis_db.Rda") # Disease DB 2019
load(file = "./gl_db.Rda") # OGlcNAc DB 2019
load(file = "./me_db.Rda") # Methylation DB 2019
load(file = "./ni_db.Rda") # Nitration DB 2019
load(file = "./p_db.Rda") # Phosphorylation DB 2019
=======
path1 <- "/Users/JCA/ptm_outdropbox/ptm/sysdata"
path2 <- "/Users/JCA/ptm_outdropbox/ptm/Rptm"

setwd(path1)

# The files corresponding to the commented lines have
# been move to GitHub.

load(file = "./aai.Rda") # Amino acids properties
# load(file = "./ac_db.Rda") # Acetylation DB 2019
# load(file = "./d.Rda") # PDB-UniProt 2019
# load(file = "./dis_db.Rda") # Disease DB 2019
# load(file = "./gl_db.Rda") # OGlcNAc DB 2019
# load(file = "./me_db.Rda") # Methylation DB 2019
# load(file = "./ni_db.Rda") # Nitration DB 2019
# load(file = "./p_db.Rda") # Phosphorylation DB 2019
>>>>>>> 6c8cdb2e66f4cc7107ac506cf5cb3c21b5df9a45
load(file = "./pax.ath.Rda") # Protein Abundance DB (Arabidopsis thaliana) 2019
load(file = "./pax.bta.Rda") # Protein Abundance DB (Bos taurus) 2019
load(file = "./pax.ddi.Rda") # Protein Abundance DB (Dictyostelium discoideum) 2019
load(file = "./pax.dme.Rda") # Protein Abundance DB (Drosophila melanogaster) 2019
load(file = "./pax.eca.Rda") # Protein Abundance DB (Equus caballus) 2019
load(file = "./pax.eco.Rda") # Protein Abundance DB (Escherichia coli) 2019
load(file = "./pax.gga.Rda") # Protein Abundance DB (Gallus gallus) 2019
load(file = "./pax.hela.Rda") # Protein Abundance DB (HeLa cells) 2019
load(file = "./pax.hsa.Rda") # Protein Abundance DB (Homo sapiens) 2019
load(file = "./pax.jurkat.Rda") # Protein Abundance DB (Jurkat cells) 2019
load(file = "./pax.mmu.Rda") # Protein Abundance DB (Mus musculus) 2019
load(file = "./pax.mtu.Rda") # Protein Abundance DB (Mycobacterium tuberculosis) 2019
load(file = "./pax.rno.Rda") # Protein Abundance DB (Rattus norvegicus) 2019
load(file = "./pax.sce.Rda") # Protein Abundance DB (Saccharomyces cerevisiae) 2019
<<<<<<< HEAD
load(file = "./reg_db.Rda") # Regulatory DB 2019
load(file = "./Sni_db.Rda") # S-nitrosylation DB 2019
load(file = "./su_db.Rda") # Sumoylation DB 2019
load(file = "./ub_db.Rda") # Ubiquitination DB 2019

setwd(path2)

usethis::use_data(aai, ac_db, d, dis_db, gl_db, me_db, ni_db, p_db, pax.ath, pax.bta, 
                  pax.ddi, pax.dme, pax.eca, pax.eco, pax.gga, pax.hela, pax.hsa, 
                  pax.jurkat, pax.mmu, pax.mtu, pax.rno, pax.sce, reg_db, Sni_db, 
                  su_db, ub_db, internal = TRUE)
=======
# load(file = "./reg_db.Rda") # Regulatory DB 2019
# load(file = "./sni_db.Rda") # S-nitrosylation DB 2019
# load(file = "./su_db.Rda") # Sumoylation DB 2019
# load(file = "./ub_db.Rda") # Ubiquitination DB 2019

setwd(path2)

usethis::use_data(aai, pax.ath, pax.bta, pax.ddi, pax.dme,
                  pax.eca, pax.eco, pax.gga, pax.hela, pax.hsa,
                  pax.jurkat, pax.mmu, pax.mtu, pax.rno, pax.sce,
                  internal = TRUE)
>>>>>>> 6c8cdb2e66f4cc7107ac506cf5cb3c21b5df9a45

# setwd(path1)
# save(aai, file = "./aai.Rda") # Amino acids properties
# save(ac_db, file = "./ac_db.Rda") # Acetylation DB 2019
# save(d, file = "./d.Rda") # PDB-UniProt 2019
# save(dis_db, file = "./dis_db.Rda") # Disease DB 2019
# save(gl_db, file = "./gl_db.Rda") # OGlcNAc DB 2019
# save(me_db, file = "./me_db.Rda") # Methylation DB 2019
# save(ni_db, file = "./ni_db.Rda") # Nitration DB 2019
# save(p_db, file = "./p_db.Rda") # Phosphorylation DB 2019
# save(pax.ath, file = "./pax.ath.Rda") # Protein Abundance DB (Arabidopsis thaliana) 2019
# save(pax.bta, file = "./pax.bta.Rda") # Protein Abundance DB (Bos taurus) 2019
# save(pax.ddi, file = "./pax.ddi.Rda") # Protein Abundance DB (Dictyostelium discoideum) 2019
# save(pax.dme, file = "./pax.dme.Rda") # Protein Abundance DB (Drosophila melanogaster) 2019
# save(pax.eca, file = "./pax.eca.Rda") # Protein Abundance DB (Equus caballus) 2019
# save(pax.eco, file = "./pax.eco.Rda") # Protein Abundance DB (Escherichia coli) 2019
# save(pax.gga, file = "./pax.gga.Rda") # Protein Abundance DB (Gallus gallus) 2019
# save(pax.hela, file = "./pax.hela.Rda") # Protein Abundance DB (HeLa cells) 2019
# save(pax.hsa, file = "./pax.hsa.Rda") # Protein Abundance DB (Homo sapiens) 2019
# save(pax.jurkat, file = "./pax.jurkat.Rda") # Protein Abundance DB (Jurkat cells) 2019
# save(pax.mmu, file = "./pax.mmu.Rda") # Protein Abundance DB (Mus musculus) 2019
# save(pax.mtu, file = "./pax.mtu.Rda") # Protein Abundance DB (Mycobacterium tuberculosis) 2019
# save(pax.rno, file = "./pax.rno.Rda") # Protein Abundance DB (Rattus norvegicus) 2019
# save(pax.sce, file = "./pax.sce.Rda") # Protein Abundance DB (Saccharomyces cerevisiae) 2019
# save(reg_db, file = "./reg_db.Rda") # Regulatory DB 2019
# save(Sni_db, file = "./Sni_db.Rda") # S-nitrosylation DB 2019
# save(su_db, file = "./su_db.Rda") # Sumoylation DB 2019
# save(ub_db, file = "./ub_db.Rda") # Ubiquitination DB 2019
