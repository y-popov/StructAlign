#include "pdb.h"
#include <math.h>

int main  (int argc, char **argv)
{
  if (argc != 12)
  {
    printf("\nUsage: %s <repr-pdb> <repr_chain> <repr-dna-chain1> <maxM> <repr-dna-chain2> <align-pdb> <align-index> <outfile> <target-prot> <target-dna1> <target_dna2>", argv[0]);
    printf("\nExample: %s 3hdd.pdb A B 20 C superpos_3hdd_1puf.pdb 1 multy.pdb F D E\n\n", argv[0]);
    exit (1);
  }
  
  char *infile1, *infile2, *maxM, *outfile;
  char chain1, chain2, dna_chain1, dna_chain2, repr_dna_chain1, repr_dna_chain2, target_prot, target_dna1, target_dna2;
  unsigned int index;

  infile1 = (char *)malloc( sizeof(char)*(strlen(argv[1])+1) );
  sscanf(argv[1],"%s", infile1);
  
  sscanf(argv[2], "%c", &chain1);
  
  maxM = (char *)malloc( sizeof(char)*(strlen(argv[4])+1) );
  sscanf(argv[4],"%s", maxM);
  
  sscanf(argv[3], "%c", &repr_dna_chain1);
  sscanf(argv[5], "%c", &repr_dna_chain2);

  infile2 = (char *)malloc( sizeof(char)*(strlen(argv[6])+1) );
  sscanf(argv[6],"%s", infile2);

  sscanf(argv[7], "%u", &index);
  if (index == 0)
  {
  	chain2 = 'E';
  	dna_chain1 = 'A';
  	dna_chain2 = 'B';
  }
  if (index == 1)
  {
  	chain2 = 'F';
  	dna_chain1 = 'C';
  	dna_chain2 = 'D';
  }
  
  outfile = (char *)malloc( sizeof(char)*(strlen(argv[8])+1) );
  sscanf(argv[8],"%s", outfile);
  
  sscanf(argv[9], "%c", &target_prot);
  sscanf(argv[10], "%c", &target_dna1);
  sscanf(argv[11], "%c", &target_dna2);
  
  
  struct atom *all_atoms_dna1 = NULL, *all_atoms_dna2 = NULL;
  struct atom *atoms_dna1 = NULL, *dna1_chain1 = NULL, *dna1_chain2 = NULL; 
  struct atom *atoms_dna2 = NULL, *dna2_chain1, *dna2_chain2;
  struct atom *atoms_prot1 = NULL, *atoms_prot2 = NULL;
  struct atom *atoms_wat1 = NULL, *atoms_wat2 = NULL;
  struct atomname *list1, *list2;
  
  char *dna_chains1 = (char *)malloc( sizeof(char)*(2+1) );
  char *prot_chains1 = (char *)malloc( sizeof(char)*(2+1) );
  char *dna_chains2 = (char *)malloc( sizeof(char)*(2+1) );
  char *prot_chains2 = (char *)malloc( sizeof(char)*(2+1) );
  
  unsigned int all_m1=0, n1=0, w1=0, i, j, all_m2=0, n2=0, w2=0;
  unsigned int m1, dna1_chain1_n=0, dna1_chain2_n=0;
  unsigned int m2, dna2_chain1_n=0, dna2_chain2_n=0;
  unsigned int maxnumber = 256;
  FILE *flow_out;
  
  //puts("Reading 1st PDB file...");
  readerPDB(infile1, &all_m1, &n1, &w1, maxnumber, &all_atoms_dna1, &dna_chains1, &atoms_prot1, &prot_chains1, &atoms_wat1, &list1); 
  //printf("...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", all_m1, n1, w1);
  if (all_m1 == 0)
    {printf("Error\nFirst structure has no DNA!");
    exit(1); }
  SelectChain(atoms_prot1, n1, &atoms_prot1, &n1, chain1);
  //printf("Atoms in selected chain: %u\n\n", n1); 
  SelectChain(all_atoms_dna1, all_m1, &dna1_chain1, &dna1_chain1_n, repr_dna_chain1);
  SelectChain(all_atoms_dna1, all_m1, &dna1_chain2, &dna1_chain2_n, repr_dna_chain2);
  atomlistmerge(&atoms_dna1, &m1, dna1_chain1, dna1_chain1_n, dna1_chain2, dna1_chain2_n);
  
  unsigned int *list_P1, *list_C11, *list_OP11, *list_OP21;
  unsigned int n_P1, n_C11, n_OP11, n_OP21;
  getAtomsNumbers(atoms_dna1, m1, &list_P1, &n_P1, "P"); 
  getAtomsNumbers(atoms_dna1, m1, &list_C11, &n_C11, "C1'");
  getAtomsNumbers(atoms_dna1, m1, &list_OP11, &n_OP11, "OP1");
  getAtomsNumbers(atoms_dna1, m1, &list_OP21, &n_OP21, "OP2");
  correctC1_P(atoms_dna1, &list_C11, &n_C11, list_P1, &n_P1);
  
  unsigned int i_max = 0;
  for (i=1; i<=n_P1; i++)
  	if (strcmp(atoms_dna1[list_P1[i]].ResNumber, maxM) == 0 && atoms_dna1[list_P1[i]].Chain == repr_dna_chain1)
  	{
  		i_max = i;
		//printf("imax:%u maxM:%s res:%s\n", i_max, maxM, atoms_dna1[list_P1[i]].ResNumber);
  	}
  if (i_max==0)
  {
  	printf("Cannot find maxM-nucleotide in %s\n", infile2);
  	exit(1);
  }

  //puts("Reading 2nd PDB file...");
  readerPDB(infile2, &all_m2, &n2, &w2, maxnumber, &all_atoms_dna2, &dna_chains2, &atoms_prot2, &prot_chains2, &atoms_wat2, &list2); 
  //printf("...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", all_m2, n2, w2);
  if (all_m2 == 0)
    {printf("Error\nSecond structure has no DNA!");
    exit(1); }
  SelectChain(atoms_prot2, n2, &atoms_prot2, &n2, chain2);
  //printf("Atoms in selected chain: %u\n\n", n2);
  SelectChain(all_atoms_dna2, all_m2, &dna2_chain1, &dna2_chain1_n, dna_chain1);
  SelectChain(all_atoms_dna2, all_m2, &dna2_chain2, &dna2_chain2_n, dna_chain2);
  atomlistmerge(&atoms_dna2, &m2, dna2_chain1, dna2_chain1_n, dna2_chain2, dna2_chain2_n);
  
  struct coordsystem dnares;
  dnares = ChangeSystem(dna1_chain1, dna1_chain1_n, &dna1_chain1, atoms_dna1[list_P1[i_max]], atoms_dna1[list_C11[i_max+1]], atoms_dna1[list_OP11[i_max]], atoms_dna1[list_OP21[i_max]], 'X'); 
  
  ChangeSystemR(atoms_prot2, n2, &atoms_prot2,  atoms_dna1[list_P1[i_max]], target_prot, dnares);  
  ChangeSystemR(dna2_chain1, dna2_chain1_n, &dna2_chain1,  atoms_dna1[list_P1[i_max]], target_dna1, dnares); 
  ChangeSystemR(dna2_chain2, dna2_chain2_n, &dna2_chain2,  atoms_dna1[list_P1[i_max]], target_dna2, dnares); 
  
  
  writetoPDB(outfile, dna2_chain1, dna2_chain1_n);
  writetoPDB(outfile, dna2_chain2, dna2_chain2_n);
  writetoPDB(outfile, atoms_prot2, n2);
  
}
