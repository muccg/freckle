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
typedef unsigned int TupleID;

/* the type of var our C and D tables are */
typedef unsigned int TupleStore;

extern const char *Bases;
extern const char *Aminos;
extern const char *TranslateUniversal;

// this is a bit mask for the tuple id
#define BASE_MASK(basestring)	(basestring==Bases?7:(basestring==Aminos?31:0))

// how many basepairs we have
#define BASE_PAIRS(basestring)	(strlen(basestring))

// this is how far to shift the bits to squeeze another in
#define BASE_BIT_SHIFT(basestring)	(basestring==Bases?3:(basestring==Aminos?5:0))		



// Function prototypes
int ipow(int x, int n);
TupleID getTupleID(const char *tuple, int len, const char *bases=Bases);
TupleStore **buildMappingTables( const char *sequence, int ktuplesize, const char *bases=Bases );
void freeMappingTables(int **tables);
int sum(int *buffer, int length);
int matchAboveThreshold(const char *seq1, int p1, const char *seq2, int p2, int k, int threshold, int window);
DotStore *doComparison(unsigned int **tables, const char *tablesequence, const char *newsequence, int ktuplesize, int window, int mismatch, int minmatch, const char *bases=Bases );
DotStore *makeDotComparison(const char *seq1, const char *seq2, int ktuplesize, int window, int mismatch, int minmatch);

// lbdot comparison
char *strrev( char *str);
void Init_code_tables();
void EncodeNTSeq(const char *seq, int p1, int p2, int *c,int *d, int nm, int nMaxDNAKtup);
int EncodeNTSeqConditional(const char *seq, int p1, int p2, int *c,int *d,int *cd, int nm, int maxHints, int nMaxDNAKtup);
int GetNtCode(const char *seq, int ktup, int intval, const int *v);
void ComplementSeq(char  *a);
char *RCseq(char *a);
DotStore **DoFastComparison(char *Seq1, char *Seq2, int SeqLen1, int SeqLen2,
						   int CompWind,int CompMism, int nMaxRepeatKtup, int nMaxDNAKtup);

// helper functions
DotStore *NewDotStore();
void DelDotStore(DotStore *store);
void DotStoreAddDot(DotStore *store, int x, int y, int len);
int DotStoreGetDotX(DotStore *store, int index);
int DotStoreGetDotY(DotStore *store, int index);
int DotStoreGetDotLength(DotStore *store, int index);
Dot *DotStoreGetDot(DotStore *store, int index);
int DotStoreGetNumDots(DotStore *store);
void DotStoreCreateIndex(DotStore *store);
void DotStoreDestroyIndex(DotStore *store);
int *DotStoreToBuffer(DotStore *store);
void DotStoreFromBuffer(DotStore *store, int *buffer);
int DotStoreBufferSize(DotStore *store, int *buffer);
void FreeIntBuffer(int *buffer);
DotStore *DotStoreFilter(DotStore *store, int minlen);
void DotStoreInterpolate(DotStore *store, int window);

// maximums
void DotStoreSetMaxX(DotStore *store, int max);
void DotStoreSetMaxY(DotStore *store, int max);
int DotStoreGetMaxX(DotStore *store);
int DotStoreGetMaxY(DotStore *store);

// conservation helper functions
Dot *DotStoreGetIndexLongestMatchingRowDot(DotStore *store, int x);
Dot *DotStoreGetIndexLongestMatchingColumnDot(DotStore *store, int y);

// DotGrid helper functions and wrappers
DotGrid *NewDotGrid();
void DelDotGrid(DotGrid *grid);
int DotGridWidth(DotGrid *grid);
int DotGridHeight(DotGrid *grid);
int DotGridGetSize(DotGrid *grid);
int DotGridGetMax(DotGrid *grid);
int DotGridGetMin(DotGrid *grid);
unsigned char *DotGridToString(DotGrid *grid);
void FreeString(unsigned char *string);
void DotGridCalculate(DotGrid *grid, DotStore *source, double x1, double y1, double x2, double y2, double scale, int window);
void DotGridAddInplace(DotGrid *source, DotGrid *add);
void DotGridFlipInplace(DotGrid *grid);
int *DotGridCalculateHistogram(DotGrid *grid);



// test debug
void getInfo();





}

#endif
