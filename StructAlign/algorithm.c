#include "pdb.h"
#include <math.h>

int main  (int argc, char **argv)
{

  /* Checking command line, 
  throw exception if not complete */

  if (argc != 4)
  {
    printf("\nUsage: %s <input file 1.pdb> <input file 2.pdb> <output file>", argv[0]);
    printf("\nExample: %s 3hdd.pdb 1puf.pdb superpos_3hdd_1puf.pdb\n\n", argv[0]);
    exit (1);
  } /* end if */
  
  /* Done checking command line */

  /* Reading command line arguments */

  char *infile1, *infile2, *outfile;

  infile1 = (char *)malloc( sizeof(char)*(strlen(argv[1])+1) );
  sscanf(argv[1],"%s", infile1);

  infile2 = (char *)malloc( sizeof(char)*(strlen(argv[2])+1) );
  sscanf(argv[2],"%s", infile2);

  outfile = (char *)malloc( sizeof(char)*(strlen(argv[3])+1) );
  sscanf(argv[3], "%s", outfile);

  /* Done reading arguments */

  /* DECLARATION AND ASSIGNMENT 
  OF SOME VARIABLES */

  struct atom *atoms_prot1 = NULL, *atoms_prot2 = NULL;
  struct atom *atoms_dna1 = NULL, *atoms_dna2 = NULL;
  struct atom *atoms_wat1 = NULL, *atoms_wat2 = NULL;
  struct atom *atoms_prot_CA1 = NULL, *atoms_prot_CA2 = NULL;
  struct atom *atoms_prot_C1 = NULL, *atoms_prot_C2 = NULL;
  struct atom *atoms_prot_i1 = NULL, *atoms_prot_i2 = NULL;
  struct atom *atoms_prot_j1 = NULL, *atoms_prot_j2 = NULL;
  struct atom *C_atoms_prot_i1 = NULL, *C_atoms_prot_j2 = NULL;

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
  
  /* Reading PDB file */
  
  puts("Reading 1st PDB file...");
  readerPDB(infile1, &m1, &n1, &w1, maxnumber, &atoms_dna1, &dna_chains1, &atoms_prot1, &prot_chains1, &atoms_wat1, &list1); 
    // pdb.c function, read PDB and put protein, DNA and water atoms in three different arrays
  printf("...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", m1, n1, w1);
  SelectChain(atoms_prot1, n1, &atoms_prot1, &n1, prot_chains1[1]); 
    // pdb.c function, select only first chain of protein from PDB file.
  printf("Atoms in selected chain: %u\n\n", n1); 
    // print to stdout number of protein atoms in selected chain

  puts("Reading 2nd PDB file...");
  readerPDB(infile2, &m2, &n2, &w2, maxnumber, &atoms_dna2, &dna_chains2, &atoms_prot2, &prot_chains2, &atoms_wat2, &list2); 
  printf("...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", m2, n2, w2);
  SelectChain(atoms_prot2, n2, &atoms_prot2, &n2, prot_chains2[1]);
  printf("Atoms in selected chain: %u\n\n", n2);

  dna_chains1[0] = ' '; dna_chains2[0] = ' ';
  printf("DNA1 chains: %s; DNA2 chains: %s\n", dna_chains1, dna_chains2);
  /* Done reading */

  /* Make lists of P, C1', OP1 atoms 
  of dna and CA atoms of protein */
  
  unsigned int *list_P1, *list_C11, *list_OP11, *list_OP21, *list_CA1, *list_P2, *list_C12, *list_OP12, *list_OP22, *list_CA2, *list_C1, *list_C2;
  unsigned int n_P1, n_C11, n_OP11, n_OP21, n_CA1, n_P2, n_C12, n_OP12, n_OP22, n_CA2, n_C1, n_C2;

  getAtomsNumbers(atoms_dna1, m1, &list_P1, &n_P1, "P"); 
    // pdb.c function, get only indexes of atoms
  getAtomsNumbers(atoms_dna1, m1, &list_C11, &n_C11, "C1'");
  getAtomsNumbers(atoms_dna1, m1, &list_OP11, &n_OP11, "OP1");
  getAtomsNumbers(atoms_dna1, m1, &list_OP21, &n_OP21, "OP2");
  getAtomsNumbers(atoms_prot1, n1, &list_CA1, &n_CA1, "CA");
  getAtomsNumbers(atoms_prot1, n1, &list_C1, &n_C1, "C");
  correctC1_P(atoms_dna1, &list_C11, &n_C11, list_P1, n_P1);
    //corrects list_C1 to use rule list_P[i] and list_C1[i+1] are in the same nucleotide

  getAtomsNumbers(atoms_dna2, m2, &list_P2, &n_P2, "P");
  getAtomsNumbers(atoms_dna2, m2, &list_C12, &n_C12, "C1'");
  getAtomsNumbers(atoms_dna2, m2, &list_OP12, &n_OP12, "OP1");
  getAtomsNumbers(atoms_dna2, m2, &list_OP22, &n_OP22, "OP2");
  getAtomsNumbers(atoms_prot2, n2, &list_CA2, &n_CA2, "CA");
  getAtomsNumbers(atoms_prot2, n2, &list_C2, &n_C2, "C");
  correctC1_P(atoms_dna2, &list_C12, &n_C12, list_P2, n_P2);

  atoms_prot_CA1 = (struct atom *)malloc( sizeof(struct atom)*(n_CA1+1) );
  for (i=1; i<=n_CA1; i++) {
    atomcpy(&atoms_prot_CA1[i], atoms_prot1[list_CA1[i]]); 
      //pdb.c function. Copy all properties of atom
	}
  atoms_prot_C1 = (struct atom *)malloc( sizeof(struct atom)*(n_C1+1) );
  for (i=1; i<=n_C1; i++) {
  	atomcpy(&atoms_prot_C1[i], atoms_prot1[list_C1[i]]);
  }
  printf("Done atoms:\tP %u\tC1 %u\tOP1 %u\tOP2 %u\tCA %u\tC %u\n",n_P1, n_C11, n_OP11, n_OP21, n_CA1, n_C1);  
  
  atoms_prot_CA2 = (struct atom *)malloc( sizeof(struct atom)*(n_CA2+1) );
  for (i=1; i<=n_CA2; i++) {
    atomcpy(&atoms_prot_CA2[i], atoms_prot2[list_CA2[i]]);
	}
  atoms_prot_C2 = (struct atom *)malloc( sizeof(struct atom)*(n_C2+1) );
  for (i=1; i<=n_C2; i++) {
  	atomcpy(&atoms_prot_C2[i], atoms_prot2[list_C2[i]]);
  }
  printf("Done atoms:\tP %u\tC1 %u\tOP1 %u\tOP2 %u\tCA %u\tC %u\n",n_P2, n_C12, n_OP12, n_OP22, n_CA2, n_C2);  

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
      ChangeSystem(atoms_prot_C1, n_C1, &C_atoms_prot_i1, atoms_dna1[list_P1[i]], atoms_dna1[list_C11[i+1]], atoms_dna1[list_OP11[i]], atoms_dna1[list_OP21[i]], 'E');
      ChangeSystem(atoms_prot_CA2, n_CA2, &atoms_prot_j2, atoms_dna2[list_P2[j]], atoms_dna2[list_C12[j+1]], atoms_dna2[list_OP12[j]], atoms_dna2[list_OP22[j]], 'F');
      ChangeSystem(atoms_prot_C2, n_C2, &C_atoms_prot_j2, atoms_dna2[list_P2[j]], atoms_dna2[list_C12[j+1]], atoms_dna2[list_OP12[j]], atoms_dna2[list_OP22[j]], 'F');
	    BidirectionalHit(atoms_prot_i1, C_atoms_prot_i1, n_CA1, atoms_prot_j2, C_atoms_prot_j2, n_CA2, &list_hit, &n_hit); 
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
  
  find_compl(atoms_dna1, list_P1, list_C11, list_OP11, list_OP21, atoms_dna_P1, n_P1, &compl1, &n_first_chain);
  find_compl(atoms_dna2, list_P2, list_C12, list_OP12, list_OP22, atoms_dna_P2, n_P2, &compl2, &m_first_chain);
  BestDiag(list_measure, n_P1, n_P2, &S_max, &i_max, &j_max, &i_start, &j_start, &i_max_measure, &j_max_measure,
  	   atoms_dna1, list_P1, atoms_dna2, list_P2, compl1, compl2, n_first_chain, m_first_chain);
  // pdb.c function 

  /* Done diagonal search */


  struct atom *atoms_dna_i1 = NULL;
  struct atom *atoms_dna_j2 = NULL;
  /* Change the system to i and j coordinates, write the alignment to file */

  ChangeSystem(atoms_prot1, n1, &atoms_prot_i1, atoms_dna1[list_P1[i_max_measure]], atoms_dna1[list_C11[i_max_measure+1]], atoms_dna1[list_OP11[i_max_measure]], atoms_dna1[list_OP21[i_max_measure]], 'E'); 
    // Note that names of chains will be changed to 'E' and 'F' by default! Modify if required.
  ChangeSystem(atoms_prot2, n2, &atoms_prot_j2, atoms_dna2[list_P2[j_max_measure]], atoms_dna2[list_C12[j_max_measure+1]], atoms_dna2[list_OP12[j_max_measure]], atoms_dna2[list_OP22[j_max_measure]], 'F');
  ChangeSystem(atoms_dna1, m1, &atoms_dna_i1, atoms_dna1[list_P1[i_max_measure]], atoms_dna1[list_C11[i_max_measure+1]], atoms_dna1[list_OP11[i_max_measure]], atoms_dna1[list_OP21[i_max_measure]], 'C'); 
  ChangeSystem(atoms_dna2, m2, &atoms_dna_j2, atoms_dna2[list_P2[j_max_measure]], atoms_dna2[list_C12[j_max_measure+1]], atoms_dna2[list_OP12[j_max_measure]], atoms_dna2[list_OP22[j_max_measure]], 'D');
    
  puts("\nFix 'nan' in pdb with this data:");
  createPDB(outfile, outfile); 
    // pdb.c function
  writetoPDB(outfile, atoms_dna_i1, m1);
  writetoPDB(outfile, 
                atoms_prot_i1, n1); 
    // pdb.c function. Write protein atoms to outfile

  writetoPDB(outfile, 
                atoms_prot_j2, n2);
  writetoPDB(outfile, atoms_dna_j2, m2);
  endPDB(outfile); 
  
  return 0;
} 
/* End main */
