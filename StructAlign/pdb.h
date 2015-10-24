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
