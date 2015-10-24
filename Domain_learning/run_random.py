import random

ETS = "A_1awc.pdb  A_1k79.pdb  A_2stt.pdb  B_1k78.pdb  C_1bc7.pdb  C_1dux.pdb  D_1k79.pdb  E_1pue.pdb  F_1k78.pdb  G_1hbx.pdb A_1k6o.pdb  A_1k7a.pdb  A_2stw.pdb  B_1mdm.pdb  C_1bc8.pdb  C_1yo5.pdb  D_1k7a.pdb  F_1dux.pdb  F_1pue.pdb  H_1hbx.pdb".split()

Homeobox = "A_1akh.pdb  A_1cqt.pdb  A_1ig7.pdb  A_1o4x.pdb  A_3hdd.pdb  B_1b72.pdb  B_1fjl.pdb  B_1le8.pdb  B_2r5y.pdb  C_1apl.pdb  C_1k61.pdb  D_1hdd.pdb  P_1lfu.pdb A_1au7.pdb  A_1du0.pdb  A_1jgg.pdb  A_1puf.pdb  A_9ant.pdb  B_1b8i.pdb  B_1hf0.pdb  B_1puf.pdb  B_2r5z.pdb  C_1e3o.pdb  C_1mnm.pdb  D_1k61.pdb  P_1nk2.pdb A_1b72.pdb  A_1fjl.pdb  A_1k61.pdb  A_1yrn.pdb  B_1akh.pdb  B_1cqt.pdb  B_1jgg.pdb  B_1yrn.pdb  B_3hdd.pdb  C_1gt0.pdb  C_1oct.pdb  D_1mnm.pdb  P_1nk3.pdb A_1b8i.pdb  A_1hf0.pdb  A_1le8.pdb  A_2hdd.pdb  B_1au7.pdb  B_1du0.pdb  B_1k61.pdb  B_2hdd.pdb  B_9ant.pdb  C_1hdd.pdb  D_1apl.pdb  P_1ahd.pdb  P_1zq3.pdb".split()

HTH3 = "3_1lmb.pdb  4_1lmb.pdb  A_1lli.pdb  A_1rio.pdb  A_6cro.pdb  B_1lli.pdb  B_1rio.pdb  L_1per.pdb  L_1rpe.pdb  L_3cro.pdb  R_1per.pdb  R_1rpe.pdb  R_3cro.pdb".split()

Fams_labels = ['ETS', 'Homeobox', 'HTH3']
Fams = [ETS, Homeobox, HTH3]

f = open('run_random.sh', 'a')

for i in range(1024):
	Fam1 = random.choice(Fams_labels)
	temp_Fams = ['ETS', 'Homeobox', 'HTH3']
	temp_Fams.remove(Fam1)
	Fam2 = random.choice(temp_Fams)

	pdb1 = random.choice( Fams[ Fams_labels.index(Fam1) ] )
	pdb2 = random.choice( Fams[ Fams_labels.index(Fam2) ] )

	f.write( "./align ../restrict_domains/{0}_new/{1} ../restrict_domains/{2}_new/{3} all_random/{4}_{3} random\n".format(Fam1, pdb1, Fam2, pdb2, pdb1.replace('.pdb', '')) )

f.write('date +%T')
f.close()
