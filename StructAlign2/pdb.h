#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

struct geompoint {
	double X;
	double Y;
	double Z;
};

struct geomvector {
	double X;
	double Y;
	double Z;
};

struct coordsystem {
	struct geompoint origin;
	struct geomvector ort_X;
	struct geomvector ort_Y;
	struct geomvector ort_Z;
};

struct atom {
   char *PDBCode;
   char *ModelNumber;
   char *AtomNumber;
   char *AtomName;
   char Alter;
   char *ResType;
   char Chain;
   char *ResNumber;
   char Insertion;
   double XCoord;
   double YCoord;
   double ZCoord;
   double Occup;
   double B;
};

struct atomname {
     char *res;
     char *id;
};

/* "readatom" reads the next atom from the (already open) PDB file;
* if key = 0 - any atom, if key = 1 - any heavy (non-hydrogen) atom, 
* key = 2 - protein heavy atom, key = 3 - DNA/RNA heavy atom */
unsigned int readatom (FILE *flowpdb, struct atom *current, 
                       int *changemodel, int key);

/* "ifprot" returns TRUE if the string is a name of an aminoacid residue
*  and FALSE otherwise */
unsigned int ifprot(char *res);

/* "ifdna" returns TRUE if the string is a name of a nucleotide
*  and FALSE otherwise */
unsigned int ifdna(char *res);

/* Copying one atom to another */
void atomcpy(struct atom *target, struct atom source);

/* Coordinates of an atom as "geompoint" structure */
struct geompoint location(struct atom a);

/* Scalar product of two vectors */
double scalar (struct geomvector v, struct geomvector u);

/* Square distance between a point A and the straight line through a point B with
*  the direction V */
double sqdistance(struct geompoint A, struct geompoint B, struct geomvector V);

/* Multiplying a vector by a number */
struct geomvector num_x_vec(double coef, struct geomvector V);

/* Shifting a point by a vector */
struct geompoint pointplusvec(struct geompoint A, struct geomvector V);

/* Vector from one point to another */
struct geomvector fromto(struct geompoint A, struct geompoint B);

/* The point on the straight line (A,V) that is
*  nearest to the line (B,W) */
struct geompoint nearest(struct geompoint A,struct geomvector V,
                         struct geompoint B,struct geomvector W);

/* Sum of two vectors */
struct geomvector sumvec(struct geomvector V, struct geomvector W);

/* Vector product of two vectors */
struct geomvector vecproduct(struct geomvector V, struct geomvector W);

/* Projection of the point B onto the line (A,V) */
struct geompoint proj(struct geompoint B, 
                      struct geompoint A, struct geomvector V);

/* Derived point = t*A + (1-t)*B */
struct geompoint derpoint(double t, struct geompoint A, struct geompoint B);

/* Vector of unit length from a non-zero vector */
struct geomvector ort(struct geomvector V);

/* Turn a vector around an axis for an angle "turn" */
struct geomvector turnvec(double turn, struct geomvector axis, 
                          struct geomvector V);

/* Changing coordinates of a vector into a new coordinate system (flag = 1) 
* or from it (flag = -1) */
struct geomvector changesystem(struct geomvector V, struct coordsystem new, 
                               int flag);
/* Reads DNA, protein and water atoms from pdb-file */
unsigned int readerPDB (char *filename, unsigned int *dnanum, unsigned int *protnum, 
              unsigned int *watnum,
              unsigned int maxnumber, struct atom **patoms_dna, char **chains_dna,
              struct atom **patoms_prot, char **chains_prot, struct atom **patoms_wat,
              struct atomname **plist);

/* Check if the char is in the string */
int inArray(char c, char *string);

/* Return the nucleotide sequence of DNA atoms*/
void Seq(struct atom *atoms, unsigned int n, char **seq, char **num_seq, unsigned int *m,  FILE *max_score);

/* Returns the array of atoms of one chain*/
unsigned int SelectChain(struct atom * atoms_from, unsigned int n_from, struct atom ** atoms_to, unsigned int * n_to, char chain);

/* Returns the number of atoms from the array of one atom type*/
unsigned int getAtomsNumbers(struct atom * atoms_from, unsigned int n, unsigned int ** list_to, unsigned int * n_to, char * atoms_type);

/* If there is different number of C1' and P atoms, it deletes redundant atoms */
unsigned int correctC1_P(struct atom *atoms_from, unsigned int **list_C, unsigned int *n_C, unsigned int *list_P, unsigned int *n_P);

/* Moves coordinate system to desired atom */
struct coordsystem ChangeSystem(struct atom * atoms_from, 
                          unsigned int n, struct atom ** atoms_to, 
						  struct atom Op, struct atom Xp, struct atom Yp, struct atom Yp2,
						  char chainname);

/* Makes comparison between Calpha atoms */
unsigned int BidirectionalHit( struct atom * atoms_i, struct atom *C_atoms_i, unsigned int n_i, struct atom * atoms_j, struct atom *C_atoms_j, unsigned int n_j, unsigned int *** list, unsigned int * n_hit);

/* M = max( 0; 4.5-rho(atom1, atom2) ) */
double Measure2_p(double *measure, unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2);

/* Tries to find a complement nucleotide */
void find_compl(struct atom *atoms1, unsigned int *list_P, unsigned int *list_C1, unsigned int *list_OP1, unsigned int *list_OP2, struct atom *atoms1_P, unsigned int n_P, unsigned int *compl, unsigned int *first_chain_length, FILE *max_score);

/* Returns a fragment with max sum */
void BestDiag(double **measures, unsigned int n, unsigned int m, double *S_max,
		unsigned int *i_max, unsigned int *j_max, unsigned int *i_start,
		unsigned int *j_start, unsigned int *i_max_measure, unsigned int *j_max_measure,
		struct atom *atoms1, unsigned int *list_P1, struct atom *atoms2, unsigned int *list_P2,
		unsigned int compl1, unsigned int compl2, unsigned int n_1chain, unsigned int m_1chain);

/* Makes a pdb-file */
unsigned int createPDB(char *filename, char *title);
unsigned int writetoPDB(char *filename, 
                struct atom *atoms, unsigned int natoms);
unsigned int endPDB(char *filename);

/* execute find_pair */
unsigned int run_3dna(char *pdb_name, unsigned int **compl, unsigned int ***compl_pairs, char ***pairs, unsigned int *len, unsigned int server, char *max_score_filename);

/* append list to another list */
void atomlistmerge(struct atom **target, unsigned int *len, struct atom *source1, unsigned int len1, struct atom *source2, unsigned int len2);

/* copies struct atom list */
void atomlistcpy(struct atom **target, struct atom *source, unsigned int leng);

/* find_compl based on 3dna data */
unsigned int run_find_compl(struct atom *atoms1_P, unsigned int n_P, unsigned int *compl, unsigned int *first_chain_length, unsigned int *compl_pair);

