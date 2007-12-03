/*
** Implementation of algorithm detailed in "Rapid and sensitive dot-matrix methods for genome analysis" by Yang and Zuang, Bioinformatics 20(4), Oxford University Press, 2004.
** Programmed by Crispin Wellington, December 2007
** Copyright Center of Comparative Genomics, Murdoch University, WA
**
** This algorithm involves breaking the DNA sequence into overlapping word pairs that are k bases long.
** These are called k-tuples. For a ktuple of size 2 and 4 base pairs we have 16 possible k-tuples.
** AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT. We say this table is l in length.
**
** We then build two tables. Table C is l units long with each array element corresponding to the relevant
** k-tuple. The values in this table at the end of the build point to the sequence position of the last 
** occurence of that k-tuple, or are 0 if that k-tuple does not appear anywhere in the sequence. While building C,
** we also build table D. Table D is just slightly shorter than the length of the sequence. If l is the
** length of the sequence then the table is l-k+1 long. Each entry in this table points (via an index) to
** the previous occurence of the relevant k-tuple whose sequence position it is in. It is set to 0 if its
** the earliest occurence of the tuple.
**
** Once these tables are complete we can take any k-tuple and quickly and efficiently find all the
** occurences in the sequence. We look at a k-tuples position in table C and that gives us the last
** occurance of that ktuple in the sequence, say at location i. Then if we look at D[i] it gives us
** the previous position j, and if we look at D[j] it gives us the previous position. In this way we
** can efficiently find every location of that k-tuple in the sequence by following these indexes
** until D[i]=0 which marks the first (we are scanning backwards) and final occurance of the tuple.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#define BASEPAIRS	4

typedef unsigned int u32;
typedef unsigned long long int u64;

extern "C" {




const char *Bases="ACGT";
#define BASE_BIT_SHIFT	2
#define BASE_MASK 3

/* this gives us a maximum ktuple size of 32 for 4 base pairs. But will probably be slower on a 32-bit system */
//typedef u64 TupleID;

/* this gives us a maximum ktuple size of 16 for 4 base pairs. But will work well on 32 bit systems */
typedef u32 TupleID;

/* find the index of the first occurence of character c in the passed in string */
// int getstringindex(const char *str, char c)
// {
// 	printf("string: %s\nchar %c\n",str,c);
// 
// 	int i=0;
// 	while(str[i] && str[i]!=c)
// 		i++;
// 	return i;
// }

/*
** ipow(x,n)
** =========
** calculates X to the power of n for integers. Uses tail recursion
*/
int ipow(int x, int n)
{
	return (!n)?1:(n&1)?x*ipow(x,n-1):ipow(x*x,n/2);
}

/*
** getTupleID
** ==========
** pass this a pointer to the base of a sequence string and a length.
** It will return the ktuple index of that sequence. It does this without building
** a lookup table to save memory.
*/
TupleID getTupleID(const char *tuple, int len)
{
	TupleID id=0;
	int index=0;

	for(int i=0; i<len; i++, tuple++)
	{
 		index=strchr(Bases, (int)*tuple)-Bases;			// get the index in Bases of where this character appears
//		index=getstringindex(Bases,*tuple);			// get the index in Bases of where this character appears
		assert(index>=0 && index<(int)strlen(Bases));		// it should be within range
		assert((index & BASE_MASK) == index);			// it shouldn't have any extra bits
		
		id=id<<BASE_BIT_SHIFT;		// right shift the id to make space for this new sequence char
	
		id |= index;			// add the index bits onto the right
	}

	// make sure it doesnt overflow
	assert(id+1);

	return id+1;				// 1 offset rather than 0
}

/*
** buildMappingTable
** =================
** Given a sequence as a string, this function builds the mapping tables C and D
**
** Input:
**   sequence:		The sequence represented as a string of characters delimited by a NULL
**   ktuplesize:	The size of the ktuple word in characters
**
** Output:
**   pointertable[2]	An array of two pointer. table[0] points to the beginning of C, and table[1] to D
*/
int **buildMappingTables( const char *sequence, int ktuplesize )
{
	int seqlen=strlen(sequence);
	assert(strlen>0);

	int ktuplearraysize=ipow(BASEPAIRS,ktuplesize);
	int darraysize=seqlen-ktuplesize+1;
	
	// we allocate our arrays
	int *C=new int [ktuplearraysize];
	int *D=new int [darraysize];

	// lets zero our arrays
	memset(C, 0, sizeof(int)*ktuplearraysize);
	memset(D, 0, sizeof(int)*darraysize);

	// Initialise D
	const char *tuple=sequence;
	for(int i=0; i<darraysize; i++, tuple++)
		//read the ith tuple from the sequence. Put its tuple id into D
		D[i]=getTupleID(tuple,ktuplesize);

	// Build the tables
	tuple=sequence;
	int cval=0, dval=0;
	for(int i=0; i<darraysize; i++, tuple++)
	{
		// the index, i, of D is assigned to C[D[i]]. So if D[1]=3, then C[3]=1
		// of course all our C array offsets are 0 based, but the algorithm REQUIRES 1 based, because 0 is a terminator
		// if C[D[i]]!=0, then it already has an index. that value should be stored in D[i] before the C value is set
		dval=D[i];
		cval=C[ dval-1 ];
		if(cval)
		{
			//has a value already
			C[dval-1]=i+1;
			D[i]=cval;
		}
		else
		{
			// is empty
			C[ dval-1 ] = i+1;
			D[i]=0;
		}
	}
	
	int **pointertable=new int *[2];
	pointertable[0]=C;
	pointertable[1]=D;

	return pointertable;
}

void freeMappingTables(int **tables)
{
	delete tables[0];
	delete tables[1];
	delete tables;
}

/*
** doComparison
** ============
** compare one sequence against the table constructed sequence
*/
void doComparison(int **tables, const char *tablesequence, const char *newsequence, int ktuplesize)
{
	
}

/*
** getInfo
** =======
** print out some info about this system
*/
void getInfo()
{
	printf("sizeof(int)=%d\n",sizeof(int));
	printf("sizeof(long int)=%d\n",sizeof(long int));
	printf("sizeof(long long int)=%d\n",sizeof(long long int));

	printf("AGC=%d\n",getTupleID("AGC",3));
	printf("GTA=%d\n",getTupleID("GTA",3));
	printf("TCA=%d\n",getTupleID("TCA",3));
	printf("AAA=%d\n",getTupleID("AAA",3));
	printf("TTT=%d\n",getTupleID("TTT",3));
	
	// from the paper itself
	//buildMappingTables("AGCTCGATCGAGTCTCGAGTAG",2);

	//buildMappingTables("GATTACAATTAACTGATCGATCGTAGCTACATGCTGACTACTGACTGCATGCATGACTGCATGCATTGACTGACTGCATGACTGCATG",8);
	int **pt=buildMappingTables("AGCTCGATCGAGTCTCGAGTAG",2);

	printf("%d - %d\n",(int)pt[0], (int)pt[1]);

	freeMappingTables(pt);
}




}

