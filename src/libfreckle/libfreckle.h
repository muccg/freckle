#ifndef _LIBFRECKLE_H_
#define _LIBFRECKLE_H_

#include "DotStore.h"
#include "DotGrid.h"

extern "C" {


typedef unsigned int u32;
typedef unsigned long long int u64;


/* this gives us a maximum ktuple size of 32 for 4 base pairs. But will probably be slower on a 32-bit system */
//typedef u64 TupleID;

/* this gives us a maximum ktuple size of 16 for 4 base pairs. But will work well on 32 bit systems */
typedef u32 TupleID;

const char *Bases="ACGT";
const char *Aminos="ACDEFGHIKLMNPQRSTVWY-";					// - = stop codon
const char TranslateUniversal[]="KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV-Y-YSSSS-CWCLFLF";

// this is a bit mask for the tuple id
#define BASE_MASK(basestring)	(basestring==Bases?3:(basestring==Aminos?31:0))

// how many basepairs we have
#define BASE_PAIRS(basestring)	(strlen(basestring))

// this is how far to shift the bits to squeeze another in
#define BASE_BIT_SHIFT(basestring)	(basestring==Bases?2:(basestring==Aminos?5:0))		



// Function prototypes
int ipow(int x, int n);
TupleID getTupleID(const char *tuple, int len, const char *bases=Bases);
int **buildMappingTables( const char *sequence, int ktuplesize, const char *bases=Bases );
void freeMappingTables(int **tables);
int sum(int *buffer, int length);
int matchAboveThreshold(const char *seq1, int p1, const char *seq2, int p2, int k, int threshold, int window);
DotStore *doComparison(int **tables, const char *tablesequence, const char *newsequence, int ktuplesize, int window, int mismatch, int minmatch, const char *bases=Bases );
DotStore *makeDotComparison(const char *seq1, const char *seq2, int ktuplesize, int window, int mismatch, int minmatch);

// helper functions
int GetDotX(DotStore *store, int index);
int GetDotY(DotStore *store, int index);
int GetDotLength(DotStore *store, int index);
Dot *GetDot(DotStore *store, int index);
int GetNumDots(DotStore *store);
void FreeDotStore(DotStore *store);
// unsigned char *DotStoreImageToString(DotStore *store, int xseqsize, int yseqsize, int longest, int window);
void DumpDotStore(DotStore *store);

// DotGrid helper functions and wrappers
DotGrid *NewDotGrid();
void DelDotGrid(DotGrid *grid);
int DotGridWidth(DotGrid *grid);
int DotGridHeight(DotGrid *grid);
int DotGridGetSize(DotGrid *grid);
int DotGridGetMax(DotGrid *grid);
int DotGridGetMin(DotGrid *grid);
unsigned char *DotGridToString(DotGrid *grid);
void DotGridCalculate(DotGrid *grid, DotStore *source, double x1, double y1, double x2, double y2, double scale, int window);
void DotGridAddInplace(DotGrid *source, DotGrid *add);
void DotGridFlipInplace(DotGrid *grid);



// test debug
void getInfo();





}

#endif
