#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pdb.h"


void read_my_pdbs(char **my_pdbs)
{
	FILE *my_pdbs_file;
	char s[9];
	int i;
	
	my_pdbs_file = fopen("my_pdbs.txt", "r");
	i = 1;
	while ( fgets(s, 10, my_pdbs_file) )
	{
		sscanf(s, "%s", s);
		strcpy( my_pdbs[i], s );
		i += 1;
	}
}

unsigned int SelectResidues(struct atom * atoms_from, unsigned int n_from, struct atom ** atoms_to, unsigned int * n_to, int from, int to)
  {
  unsigned int i, j=0;
  int Resnumber;

  (*atoms_to) = (struct atom *)malloc( sizeof(struct atom)*(n_from+1) );
  /*printf("OK1\n");*/ // Enable in test mode
  for (i=1; i<=n_from; i++) {
    sscanf(atoms_from[i].ResNumber, "%i", &Resnumber);
    if ( Resnumber >= from  && Resnumber <= to) {
	  j++;
	  (*atoms_to)[j] = atoms_from[i];
	  }
	}
  /*printf("OK\n");*/ // Enable in test mode
  *n_to = j;
    
  return 0;
}

void main()
{

FILE *domains_file;
char *pdb_file_name;
pdb_file_name = (char *)malloc( sizeof(char)*100 );
char string[125];
char *pdb;
char temp[5];
char chain;
pdb = (char *)malloc( sizeof(char)*10 );

int i, start, end, flag;
unsigned int n1, m1;
unsigned int maxnumber = 256;

struct atom *atoms_prot1 = NULL;
struct atom *atoms_dna1 = NULL;
struct atomname *list1;

char **my_pdbs;
my_pdbs = (char **)malloc( sizeof(char *)*(60+1) ); //edit number of my pdb
for (i=1; i<=60; i++)
	{ my_pdbs[i] = (char *)malloc( sizeof(char)*(8+1) ); }


domains_file = fopen("Homeobox_domains.txt", "r");
if (domains_file == NULL)
{
	puts("No such out file!");
	exit(1);
}

read_my_pdbs(my_pdbs);

while ( fgets(string, 125, domains_file) )
{
	start = 0;
	end = 0;
	if ( string[6]=='_' )
		{ sscanf(string, "%*s %s %c:", pdb, &chain); }
	else
		{ sscanf(string, "%*s %s %c:%d-%d", pdb, &chain, &start, &end); }
	
	flag = 0;
	strcat(pdb, ".pdb");
	for (i=1; i<=60; i++) //edit number of my pdb
	{
		//printf("%s -- %s\n", my_pdbs[i], pdb);
		if ( strcmp(my_pdbs[i], pdb)==0 )
			{ flag = 1; break; }
	}
	if (flag == 1)
	{
		strcpy(pdb_file_name, "../pdb/pdbs/Homeobox/");
		strcat(pdb_file_name, pdb);
		readerPDB(pdb_file_name, &m1, &n1, maxnumber, &atoms_dna1, &atoms_prot1, &list1); 
		SelectChain(atoms_prot1, n1, &atoms_prot1, &n1, chain); 
		if (start != 0 && end != 0)
		{
			SelectResidues(atoms_prot1, n1, &atoms_prot1, &n1, start, end);
		}
		strcpy(pdb_file_name, "Homeobox_new/");
		sprintf(temp, "%c_", chain);
		strcat(pdb_file_name, temp);
		strcat(pdb_file_name, pdb);
		createPDB(pdb_file_name, pdb);
		writetoPDB(pdb_file_name, atoms_dna1, m1);
		writetoPDB(pdb_file_name, atoms_prot1, n1);
		printf("File %s was created\n", pdb_file_name);
		
	}
}


}
