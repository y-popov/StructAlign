#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pdb.h"

void process_pdb(char *pdb, FILE *out_file)
{
	unsigned int d, d_chain1, d_chain2; //numbers of atoms
	unsigned int p, i, j, c;
	unsigned int maxnumber = 256;
	unsigned int dna_chain_num, compl_P;
	int k; //+1 or -1: direction for mooving
	char user_choice, compl_chain;
	char line[5], compl_P_string[5];
	char *dna_chains = (char *)malloc( sizeof(char)*(2+1) );
	struct atom *atoms_dna = NULL; //for all dnas atoms
	struct atom *atoms1 = NULL; //for atoms in the selected chain
	struct atom *atoms2 = NULL;
	struct atom *changed_P = NULL;
	struct atomname *list;

	readerPDB(pdb, &d, &p, maxnumber, &atoms_dna, &dna_chains, &dna_chain_num, &list);
	printf("%u\n", dna_chain_num);
	
	for (i=1; i<=dna_chain_num; i++)
	{
		printf("%c\n", dna_chains[i]);
		SelectChain(atoms_dna, d, &atoms1, &d_chain1, dna_chains[i]);
		printf("Atoms in the %c chain: %u\n", dna_chains[i], d_chain1);
		
		unsigned int *list_P1, *list_P2, *list_C1, *list_OP1, *list_OP2;
		unsigned int n_P1, n_P2, n_C1, n_OP1, n_OP2;
		getAtomsNumbers(atoms1, d_chain1, &list_P1, &n_P1, "P"); 
		getAtomsNumbers(atoms1, d_chain1, &list_C1, &n_C1, "C1'"); 
		getAtomsNumbers(atoms1, d_chain1, &list_OP1, &n_OP1, "OP1"); 
		getAtomsNumbers(atoms1, d_chain1, &list_OP2, &n_OP2, "OP2"); 
		correctC1_P(atoms1, &list_C1, &n_C1, list_P1, n_P1);
		
		printf("Enter a complement chain to chain %c: ", dna_chains[i]);
		if (fgets(line, sizeof(line), stdin) == NULL) 
		{
		        fprintf(stderr, "Error!\n");
		        exit(1);
		}
		sscanf(line, "%c", &compl_chain);
		compl_chain = toupper(compl_chain);
		printf("Your input: %c\n", compl_chain);
		
		c = 0; //counter for first chain
		do
		{
			c++;
			printf("First P atom is in %s nucleotide\n", atoms1[ list_P1[c] ].ResNumber);
			printf("Enter a number of a complement nucleotide: ");
			if (fgets(line, sizeof(line), stdin) == NULL) 
			{
				fprintf(stderr, "Error!\n");
				exit(1);
			}
			sscanf(line, "%u", &compl_P);
			printf("Your input: %u\n", compl_P);
		}
		while (compl_P == 0);
		
		SelectChain(atoms_dna, d, &atoms2, &d_chain2, compl_chain);
		getAtomsNumbers(atoms2, d_chain2, &list_P2, &n_P2, "P"); 
		
		j = 1; //counter for complement chain
		sprintf(compl_P_string, "%u", compl_P);
		while ( strcmp(atoms2[ list_P2[j] ].ResNumber , compl_P_string) != 0)
			{j++;}
		if (j>3)
			{k = -1;}
		else
			{k = 1;}
		while ( (j<=n_P2) && (j>=1) && (c<=n_P1))
		{
			printf("%s--%s\n", atoms1[list_P1[c]].ResNumber, atoms2[list_P2[j]].ResNumber);
			printf("I use P in %s, C1' in %s, OP1 in %s, OP2 in %s\n", atoms1[list_P1[c]].ResNumber, atoms1[list_C1[c+1]].ResNumber, atoms1[list_OP1[c]].ResNumber, atoms1[list_OP2[c]].ResNumber);
			ChangeSystem( atoms2[list_P2[j]], 1, &changed_P, atoms1[list_P1[c]], atoms1[list_C1[c+1]], atoms1[list_OP1[c]], atoms1[list_OP2[c]], 'E'); 
			fprintf( out_file, "%s%c %s%c\t%lg\t%lg\t%lg\n",  atoms1[list_P1[c]].ResNumber,  atoms1[list_P1[c]].Chain, changed_P[1].ResNumber, changed_P[1].Chain, changed_P[1].XCoord, changed_P[1].YCoord, changed_P[1].ZCoord );
			c++; 
			j += k;
		}
	}

return;	
}

void main()
{
	char *pdb_name="1KX3.pdb";
	char *sys_cmd = (char *)malloc( sizeof(char)*(15+7+1) );
	sprintf(sys_cmd, "jmol %s &", pdb_name);
	//system(sys_cmd);
	
	FILE *outfile;
	outfile = fopen("compl_P_coords2.txt", "a");
	if (outfile == NULL)
	{
		puts("No such out file!");
		exit(1);
	}
	
	process_pdb(pdb_name, outfile);
	
	fclose(outfile);
	exit(0);
}
