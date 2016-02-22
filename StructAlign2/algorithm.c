#include "pdb.h"
#include <math.h>
#define SERVER 0

int main  (int argc, char **argv)
{
  puts("Program starts");
  /* Checking command line, 
  throw exception if not complete */

  if (argc != 7)
  {
    printf("\nUsage: %s <input file 1.pdb> <input file 2.pdb> <output file> <chain1> <chain2> <file name for score>", argv[0]);
    printf("\nExample: %s 3hdd.pdb 1puf.pdb superpos_3hdd_1puf.pdb A B max_score.txt\n\n", argv[0]);
    exit (1);
  } /* end if */
  
  /* Done checking command line */

  /* Reading command line arguments */

  char *infile1, *infile2, *outfile, *max_score_filename;
  char chain1, chain2;

  infile1 = (char *)malloc( sizeof(char)*(strlen(argv[1])+1) );
  sscanf(argv[1],"%s", infile1);

  infile2 = (char *)malloc( sizeof(char)*(strlen(argv[2])+1) );
  sscanf(argv[2],"%s", infile2);

  outfile = (char *)malloc( sizeof(char)*(strlen(argv[3])+1) );
  sscanf(argv[3], "%s", outfile);

  sscanf(argv[4], "%c", &chain1);
  sscanf(argv[5], "%c", &chain2);

  max_score_filename = (char *)malloc( sizeof(char)*(strlen(argv[6])+1) );
  sscanf(argv[6], "%s", max_score_filename);

  /* Done reading arguments */

  /* DECLARATION AND ASSIGNMENT 
  OF SOME VARIABLES */

  struct atom *atoms_prot1 = NULL, *atoms_prot2 = NULL;
  struct atom *atoms_dna1 = NULL, *atoms_dna2 = NULL;
  struct atom *atoms_wat1 = NULL, *atoms_wat2 = NULL;
  struct atom *atoms_prot_CA1 = NULL, *atoms_prot_CA2 = NULL;
  struct atom *atoms_prot_i1 = NULL, *atoms_prot_i2 = NULL;
  struct atom *atoms_prot_j1 = NULL, *atoms_prot_j2 = NULL;

  struct atomname *list1, *list2;

  unsigned int maxnumber = 256; //primary number of atoms in each groups (DNA, protein, water)
  FILE *flow_out;
  unsigned int m1=0, n1=0, w1=0, i, j, m2=0, n2=0, w2=0;
  /* m - number of DNA atoms
   * n - number of protein atoms
   * w - number of water atoms */
  char *dna_chains1 = (char *)malloc( sizeof(char)*(2+1) );
  char *prot_chains1 = (char *)malloc( sizeof(char)*(2+1) );
  char *dna_chains2 = (char *)malloc( sizeof(char)*(2+1) );
  char *prot_chains2 = (char *)malloc( sizeof(char)*(2+1) );

  /*** For server ***/
  FILE *max_score;

//printf("%s", max_score_filename);
  max_score = fopen(max_score_filename, "w");
  if (max_score == NULL) {
    perror(max_score_filename);
    exit(1);
  }
  
  /* Reading PDB file */
  
  puts("Reading 1st PDB file...");
  readerPDB(infile1, &m1, &n1, &w1, maxnumber, &atoms_dna1, &dna_chains1, &atoms_prot1, &prot_chains1, &atoms_wat1, &list1); 
    // pdb.c function, read PDB and put protein, DNA and water atoms in three different arrays
  printf("...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", m1, n1, w1);
  if (m1 == 0)
    {fprintf(max_score, "Error\nFirst structure has no DNA!");
    exit(1); }
  if (chain1=='@')
    chain1 = prot_chains1[1];
  else
    if ( inArray(chain1, prot_chains1)==0 )
      { fprintf(max_score, "Error\nChain %c is not a protein chain in first structure", chain1);
      exit(1); }
  SelectChain(atoms_prot1, n1, &atoms_prot1, &n1, chain1);
  // pdb.c function
  printf("Atoms in selected chain: %u\n\n", n1); 
    // print to stdout number of protein atoms in selected chain

  puts("Reading 2nd PDB file...");
  readerPDB(infile2, &m2, &n2, &w2, maxnumber, &atoms_dna2, &dna_chains2, &atoms_prot2, &prot_chains2, &atoms_wat2, &list2); 
  printf("...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", m2, n2, w2);
  if (m2 == 0)
    {fprintf(max_score, "Error\nSecond structure has no DNA!");
    exit(1); }
  if (chain2=='@')
    chain2 = prot_chains2[1];
  else
    if ( inArray(chain2, prot_chains2)==0 )
      { fprintf(max_score, "Error\nChain %c is not a protein chain in second structure!", chain2);
      exit(1); }
  SelectChain(atoms_prot2, n2, &atoms_prot2, &n2, chain2);
  printf("Atoms in selected chain: %u\n\n", n2);

  dna_chains1[0] = ' '; dna_chains2[0] = ' ';
  /*** For server ***/
  char dna_chain1, dna_chain2;
  dna_chain2 = dna_chains2[1];
  
  printf("DNA1 chains: %s; DNA2 chains: %s\n", dna_chains1, dna_chains2);
  /* Done reading */

  /*** For server: *************
  *****reading DNA sequences ***/
  struct atom *dna1_chain1, *dna1_chain2, *dna2_chain1, *dna2_chain2;
  unsigned int dna1_chain1_n=0, dna1_chain2_n=0, dna2_chain1_n=0, dna2_chain2_n=0;
  unsigned int dna_n11, dna_n12, dna_n21, dna_n22; //sequence length
  char *dna_seq11, *dna_seq12, *dna_seq21, *dna_seq22;
  char dna1_chain1_start[5], dna1_chain1_end[5], dna1_chain2_start[5], dna1_chain2_end[5], dna2_chain1_start[5], dna2_chain1_end[5], dna2_chain2_start[5], dna2_chain2_end[5];
  	SelectChain(atoms_dna1, m1, &dna1_chain1, &dna1_chain1_n, dna_chains1[1]);
  strcpy(dna1_chain1_start, dna1_chain1[1].ResNumber);
  strcpy(dna1_chain1_end, dna1_chain1[dna1_chain1_n].ResNumber);
  Seq(dna1_chain1, dna1_chain1_n, &dna_seq11, &dna_n11, max_score);

  	SelectChain(atoms_dna1, m1, &dna1_chain2, &dna1_chain2_n, dna_chains1[2]);
  strcpy(dna1_chain2_start, dna1_chain2[1].ResNumber);
  strcpy(dna1_chain2_end, dna1_chain2[dna1_chain2_n].ResNumber);
  Seq(dna1_chain2, dna1_chain2_n, &dna_seq12, &dna_n12, max_score);

  	SelectChain(atoms_dna2, m2, &dna2_chain1, &dna2_chain1_n, dna_chains2[1]);
  strcpy(dna2_chain1_start, dna2_chain1[1].ResNumber);
  strcpy(dna2_chain1_end, dna2_chain1[dna2_chain1_n].ResNumber);
  Seq(dna2_chain1, dna2_chain1_n, &dna_seq21, &dna_n21, max_score);

  	SelectChain(atoms_dna2, m2, &dna2_chain2, &dna2_chain2_n, dna_chains2[2]);
  strcpy(dna2_chain2_start, dna2_chain2[1].ResNumber);
  strcpy(dna2_chain2_end, dna2_chain2[dna2_chain2_n].ResNumber);
  Seq(dna2_chain2, dna2_chain2_n, &dna_seq22, &dna_n22, max_score);
  

  /* Make lists of P, C1', OP1 atoms 
  of dna and CA atoms of protein */
  
  unsigned int *list_P1, *list_C11, *list_OP11, *list_OP21, *list_CA1, *list_P2, *list_C12, *list_OP12, *list_OP22, *list_CA2;
  unsigned int n_P1, n_C11, n_OP11, n_OP21, n_CA1, n_P2, n_C12, n_OP12, n_OP22, n_CA2;

  getAtomsNumbers(atoms_dna1, m1, &list_P1, &n_P1, "P"); 
    // pdb.c function, get only indexes of atoms
  getAtomsNumbers(atoms_dna1, m1, &list_C11, &n_C11, "C1'");
  getAtomsNumbers(atoms_dna1, m1, &list_OP11, &n_OP11, "OP1");
  getAtomsNumbers(atoms_dna1, m1, &list_OP21, &n_OP21, "OP2");
  getAtomsNumbers(atoms_prot1, n1, &list_CA1, &n_CA1, "CA");
  correctC1_P(atoms_dna1, &list_C11, &n_C11, list_P1, n_P1);
    //corrects list_C1 to use rule list_P[i] and list_C1[i+1] are in the same nucleotide

  getAtomsNumbers(atoms_dna2, m2, &list_P2, &n_P2, "P");
  getAtomsNumbers(atoms_dna2, m2, &list_C12, &n_C12, "C1'");
  getAtomsNumbers(atoms_dna2, m2, &list_OP12, &n_OP12, "OP1");
  getAtomsNumbers(atoms_dna2, m2, &list_OP22, &n_OP22, "OP2");
  getAtomsNumbers(atoms_prot2, n2, &list_CA2, &n_CA2, "CA");
  correctC1_P(atoms_dna2, &list_C12, &n_C12, list_P2, n_P2);

  atoms_prot_CA1 = (struct atom *)malloc( sizeof(struct atom)*(n_CA1+1) );
  for (i=1; i<=n_CA1; i++) {
    atomcpy(&atoms_prot_CA1[i], atoms_prot1[list_CA1[i]]); 
      //pdb.c function. Copy all properties of atom
	}
  printf("Done atoms:\tP %u\tC1 %u\tOP1 %u\tOP2 %u\tCA %u\n",n_P1, n_C11, n_OP11, n_OP21, n_CA1);  
  
  atoms_prot_CA2 = (struct atom *)malloc( sizeof(struct atom)*(n_CA2+1) );
  for (i=1; i<=n_CA2; i++) {
    atomcpy(&atoms_prot_CA2[i], atoms_prot2[list_CA2[i]]);
	}
  printf("Done atoms:\tP %u\tC1 %u\tOP1 %u\tOP2 %u\tCA %u\n",n_P2, n_C12, n_OP12, n_OP22, n_CA2);  

  /* Done making lists */

  /* Create array of measures for 
  all pairs of dna P atoms (list_measure). 
  dim(list_measure) = n_P2 x n_P1 */
  
  double **list_measure;
  unsigned int n_hit; 
  unsigned int **list_hit;
  
  list_measure = (double **)malloc( (n_P2+1)*sizeof(double *)*(n_P1+1) );
  for (i=1; i<=n_P1; i++)  list_measure[i] = (double *)malloc( sizeof(double)*(n_P2+1) );
  //for (i=1; i<=n_P1; i++) for (j=1; j<=n_P2;j++)  list_measure[i][j] = (unsigned int *)malloc( sizeof(unsigned int) ); 
    // Enable if measure function returns unsigned int 
  
  
  for (i=1; i<=n_P1; i++){
    for (j=1; j<=n_P2; j++){  
      ChangeSystem(atoms_prot_CA1, n_CA1, &atoms_prot_i1, atoms_dna1[list_P1[i]], atoms_dna1[list_C11[i+1]], atoms_dna1[list_OP11[i]], atoms_dna1[list_OP21[i]], 'E'); 
        // pdb.c function. Change the coordinate system of protein with given nucleotide
      ChangeSystem(atoms_prot_CA2, n_CA2, &atoms_prot_j2, atoms_dna2[list_P2[j]], atoms_dna2[list_C12[j+1]], atoms_dna2[list_OP12[j]], atoms_dna2[list_OP22[j]], 'F');
	    BidirectionalHit(atoms_prot_i1, n_CA1, atoms_prot_j2, n_CA2, &list_hit, &n_hit); 
        // pdb.c function.
	    Measure2_p(&(list_measure[i][j]), list_hit, n_hit, atoms_prot_i1, atoms_prot_j2); 
        // pdb.c function.
	    // printf("Measure: %f  i: %u j: %u n_P1: %u  n_P2: %u \n",(list_measure[i][j]), i,j, n_P1, n_P2); 
        // Enable in test mode
    }
  }

  /* Done creation of array of measures */


  /* Start working with diagonals */

  unsigned int i_max, j_max, i_start, j_start, i_max_measure, j_max_measure;
  unsigned int compl1, compl2, n_first_chain, m_first_chain;
  double S_max;
  struct atom *atoms_dna_P1 = NULL;
  struct atom *atoms_dna_P2 = NULL;

  // print measure-table
  /*printf("  ");
  for (j=1; j<=n_P2; j++) printf("%4d", j);
  puts("");
  for (i=n_P1; i>=1; i--){
  	printf("%2d", i);
	for (j=1; j<=n_P2; j++){
		printf("%4.0f", list_measure[i][j]>0 ? list_measure[i][j] : 0);
	}
	puts("");
  }
  printf("  ");
  for (j=1; j<=n_P2; j++) printf("%4d", j);
  puts(""); 
  */
 
  i_max_measure=0;
  j_max_measure=0;
  atoms_dna_P1 = (struct atom *)malloc( sizeof(struct atom)*(n_P1+1) );
	for (i=1; i<=n_P1; i++) {
	atomcpy(&atoms_dna_P1[i], atoms_dna1[list_P1[i]]);
  }
  atoms_dna_P2 = (struct atom *)malloc( sizeof(struct atom)*(n_P2+1) );
	for (i=1; i<=n_P2; i++) {
	atomcpy(&atoms_dna_P2[i], atoms_dna2[list_P2[i]]);
  }
  
  find_compl(atoms_dna1, list_P1, list_C11, list_OP11, list_OP21, atoms_dna_P1, n_P1, &compl1, &n_first_chain, max_score);
  find_compl(atoms_dna2, list_P2, list_C12, list_OP12, list_OP22, atoms_dna_P2, n_P2, &compl2, &m_first_chain, max_score);
  BestDiag(list_measure, n_P1, n_P2, &S_max, &i_max, &j_max, &i_start, &j_start, &i_max_measure, &j_max_measure,
  	   atoms_dna1, list_P1, atoms_dna2, list_P2, compl1, compl2, n_first_chain, m_first_chain);
  // pdb.c function 


  /* Done diagonal search */

  /*** For server work - DNA alignment making ***/
  
  unsigned int i_max_measure_compl, j_max_measure_compl;
  
  outfile[strlen(outfile)-11] = chain1;
  outfile[strlen(outfile)-5] = chain2;
  i_max_measure_compl = ((i_max_measure > n_first_chain) ? compl1-i_max_measure+1 : compl1-i_max_measure+1);
  j_max_measure_compl = ((j_max_measure > m_first_chain) ? compl2-j_max_measure+1 : compl2-j_max_measure+1);
  //printf("%c %c %c %c %c %c\n%lg\n%s\n%s", chain1, chain2, dna_chains1[1], dna_chains1[2], dna_chains2[1], dna_chains2[2], S_max, dna_string1, dna_string2);

  fprintf(max_score, "%c %c %c %c %c %c %s %s %s %s %s %s %s %s %s %s %s %s\n%lg\n%s\n%s\n%s\n%s", chain1, chain2, dna_chains1[1], dna_chains1[2], dna_chains2[1], dna_chains2[2], dna1_chain1_start, dna1_chain1_end, dna1_chain2_start, dna1_chain2_end, dna2_chain1_start, dna2_chain1_end, dna2_chain2_start, dna2_chain2_end, atoms_dna1[list_P1[i_max_measure]].ResNumber, atoms_dna2[list_P2[j_max_measure]].ResNumber, atoms_dna1[list_P1[i_max_measure_compl]].ResNumber, atoms_dna2[list_P2[j_max_measure_compl]].ResNumber, S_max, dna_seq11, dna_seq12, dna_seq21, dna_seq22);
  fclose(max_score);
  

  struct atom *atoms_dna_i1 = NULL;
  struct atom *atoms_dna_i2 = NULL;
  struct atom *atoms_dna_j1 = NULL;
  struct atom *atoms_dna_j2 = NULL;
  /* Change the system to i and j coordinates, write the alignment to file */

  ChangeSystem(atoms_prot1, n1, &atoms_prot_i1, atoms_dna1[list_P1[i_max_measure]], atoms_dna1[list_C11[i_max_measure+1]], atoms_dna1[list_OP11[i_max_measure]], atoms_dna1[list_OP21[i_max_measure]], 'E'); 
    // Note that names of chains will be changed to 'E' and 'F' by default! Modify if required.
  ChangeSystem(atoms_prot2, n2, &atoms_prot_j2, atoms_dna2[list_P2[j_max_measure]], atoms_dna2[list_C12[j_max_measure+1]], atoms_dna2[list_OP12[j_max_measure]], atoms_dna2[list_OP22[j_max_measure]], 'F');
  ChangeSystem(dna1_chain1, dna1_chain1_n, &atoms_dna_i1, atoms_dna1[list_P1[i_max_measure]], atoms_dna1[list_C11[i_max_measure+1]], atoms_dna1[list_OP11[i_max_measure]], atoms_dna1[list_OP21[i_max_measure]], 'A'); 
  ChangeSystem(dna1_chain2, dna1_chain2_n, &atoms_dna_i2, atoms_dna1[list_P1[i_max_measure]], atoms_dna1[list_C11[i_max_measure+1]], atoms_dna1[list_OP11[i_max_measure]], atoms_dna1[list_OP21[i_max_measure]], 'B'); 
  ChangeSystem(dna2_chain1, dna2_chain1_n, &atoms_dna_j1, atoms_dna2[list_P2[j_max_measure]], atoms_dna2[list_C12[j_max_measure+1]], atoms_dna2[list_OP12[j_max_measure]], atoms_dna2[list_OP22[j_max_measure]], 'C');
  ChangeSystem(dna2_chain2, dna2_chain2_n, &atoms_dna_j2, atoms_dna2[list_P2[j_max_measure]], atoms_dna2[list_C12[j_max_measure+1]], atoms_dna2[list_OP12[j_max_measure]], atoms_dna2[list_OP22[j_max_measure]], 'D');
    
  //puts("\nFix 'nan' in pdb with this data:");
  createPDB(outfile, outfile); 
    // pdb.c function
  writetoPDB(outfile, atoms_dna_i1, dna1_chain1_n);
  writetoPDB(outfile, atoms_dna_i2, dna1_chain2_n);
  writetoPDB(outfile, 
                atoms_prot_i1, n1); 
    // pdb.c function. Write protein atoms to outfile

  writetoPDB(outfile, 
                atoms_prot_j2, n2);
  writetoPDB(outfile, atoms_dna_j1, dna2_chain1_n);
  writetoPDB(outfile, atoms_dna_j2, dna2_chain2_n);
  endPDB(outfile); 
  
  return 0;
} 
/* End main */
