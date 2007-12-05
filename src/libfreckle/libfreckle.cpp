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

#include "libfreckle.h"


extern "C" {

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
TupleID getTupleID(const char *tuple, int len, const char *bases)
{
	TupleID id=0;
	int index=0;

	for(int i=0; i<len; i++, tuple++)
	{
		index=strchr(bases, (int)*tuple)-bases;			// get the index in Bases of where this character appears
		assert(index>=0 && index<(int)strlen(bases));		// it should be within range
		assert((index & BASE_MASK(bases)) == index);			// it shouldn't have any extra bits
		
		id=id<<BASE_BIT_SHIFT(bases);		// right shift the id to make space for this new sequence char
	
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
int **buildMappingTables( const char *sequence, int ktuplesize, const char *bases )
{
	int seqlen=strlen(sequence);
	assert(seqlen>0);

	int ktuplearraysize=ipow(BASE_PAIRS(bases),ktuplesize);
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
		D[i]=getTupleID(tuple,ktuplesize,bases);

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
** sum function just returns the sum of a buffer of ints
*/
int sum(int *buffer, int length)
{
	int sum=0;

	for(int i=0; i<length; i++)
		sum+=buffer[i];

	return sum;
}

/*
** matchAboveThreshold
** ===================
** given two sequences, seq1 and seq2, and positions in those sequences p1 and p2, compute how long the
** two sequences match for. this assumes that a length of 'k' matches already (the tuple size)
*/
int matchAboveThreshold(const char *seq1, int p1, const char *seq2, int p2, int k, int mismatch, int window)
{
	assert(window>0);
	//assert(k>0);				//k=0 for methods 2 and 3 where we are using a coded ktuple to initiate these search locations
	int ringbuf[window];
	
	memset(ringbuf,0,sizeof(ringbuf));

	int matchlength=k;

	const char *s1=seq1+p1+k;
	const char *s2=seq2+p2+k;
	
	while((sum(ringbuf,window)<=mismatch) && *s1 && *s2 )
		ringbuf[matchlength++%window] = *s1++==*s2++?0:1;		//compare s1 and s2 characters. assign 1 to the relevant 
		
	return matchlength-((sum(ringbuf,window)<=mismatch)?0:1);
		
}


/*
** doComparison
** ============
** compare one sequence against the table constructed sequence
**
*/
DotStore *doComparison(int **tables, const char *tablesequence, const char *newsequence, int ktuplesize, int window, int mismatch, int minmatch, const char *bases )
{
	assert(mismatch<window);
	assert(window>=ktuplesize);

	int *C, *D;
	DotStore *dotstore=new DotStore();
	C=tables[0];
	D=tables[1];

	int newseqlen=strlen(newsequence);
	assert(newseqlen>0);

	int darraysize=newseqlen-ktuplesize+1;

	// go through each k-tuple on the newsequence
	const char *tuple=newsequence;
	int tupleid=0;
	for(int i=0; i<darraysize; i++, tuple++)
	{
		// first we get the id of this tuple
		tupleid=getTupleID(tuple,ktuplesize,bases);
		
		//now we look it up in the table C to find the last occurance, and move backwards 
		//through the linked list expressed in table D
		for(int position=C[tupleid-1]; position; position=D[position-1])
		{
			// so position is a tuple position match in the tabled sequence
			// now we search forward to see how long the match is (with threshold)
			int matchlen=matchAboveThreshold(tablesequence,position-1,newsequence,i,ktuplesize,mismatch,window);
			if(matchlen>=minmatch)
				dotstore->AddDot(position-1,i,matchlen);
		}
	}
	
	return dotstore;
}

/*
** makeDotComparison
** =================
** do a complete comparison including building table and comparing. returns the dotstore.
*/
DotStore *makeDotComparison(const char *seq1, const char *seq2, int ktuplesize, int window, int mismatch, int minmatch)
{
	return doComparison(buildMappingTables(seq1,ktuplesize), seq1, seq2, ktuplesize, window, mismatch, minmatch);
}

/*
** translateSequence
** =================
** translate a DNA sequence into an amino acid sequence. Returns three sequences which correspond to the
** three reading frames of the original sequence
*/
char **convertSequence(const char *sequence)
{
	printf("convertSequence()\n");
	int seqlen=strlen(sequence);

	//translation table
	TupleID codetablelen=strlen(TranslateUniversal);
	assert(codetablelen==64);

	// storage area for results. 3 strings. leach of length  floor(seqlen/3). Each is zero terminated including the initial list
	char **results=new char *[4];
	memset(results,0,sizeof(char*)*4);
	for(int i=0; i<3; i++)
	{
		results[i]=new char[seqlen/3+1];
		memset(results[i],0,sizeof(char)*(seqlen/3+1));
	}

	printf("seqlen=%d, %d\n",seqlen,seqlen/3);
	for(int i=0; i<seqlen/3; i++)
		for (int offset=0; offset<3; offset++)
		{
			if(i*3+offset+3<=seqlen)
			{
				printf("i:%d offset:%d seq[n]=%s\n",i,offset, sequence+i*3);
				TupleID id=getTupleID(sequence+i*3+offset, 3)-1;
				assert(id>=0);
				assert(id<codetablelen);
				results[offset][i]=TranslateUniversal[id];
			}
			else
			{	
				printf("overhang\n");
			}
		}

	return results;
}

/*
** makeDotComparisonByTranslation
** ==============================
** do a complete comparison including translating, building table and comparing. returns the dotstore.
** this uses method 2 of the Huang Zhang paper
*/
DotStore *makeDotComparisonByTranslation(const char *seq1, const char *seq2, int ktuplesize, int window, int mismatch, int minmatch)
{
	char **seq1translated=convertSequence(seq1);
	char **seq2translated=convertSequence(seq2);

	//build three mapping tables for the 3 reading frames of sequence 1
	int ***mappingtables=new int **[3];
	for(int i=0; i<3; i++)
	{
		printf("=%d=\n",i);
		mappingtables[i]=buildMappingTables(seq1translated[0],ktuplesize,Aminos);
	}

	printf("mapping tables done\n");
	
	//now for each amino acid match (in the 3x3 reading frames) we scan the original sequence
	//for the matchlength.
	DotStore *ds=NULL;
	DotStore *dotstore=new DotStore();
	for(int xframe=0; xframe<3; xframe++)
		for(int yframe=0; yframe<3; yframe++)
		{
			printf("doComparison %d %d...\n",xframe,yframe);
			printf("s1:%s\n",seq1translated[xframe]);
			printf("s2:%s\n",seq2translated[yframe]);
			ds=doComparison( mappingtables[xframe], seq1translated[xframe], seq2translated[yframe], ktuplesize, window, mismatch/3, 1,Aminos);
			printf("done. %d matches\n",ds->GetNum());
			for(int num=0; num<ds->GetNum(); num++)
			{
				Dot *dot=ds->GetDot(num);
				int dx=dot->x;
				int dy=dot->y;
				int length=dot->length;

				// the positions in the original sequences
				int originalx=dx*3+xframe;
				int originaly=dy*3+yframe;

				// get the original match length
				int matchlen=matchAboveThreshold(seq1, originalx, seq2, originaly, 0, mismatch, window);
				if(matchlen>=minmatch)
					dotstore->AddDot(originalx,originaly,matchlen);
				
			}
			delete ds;
		}
	
	return dotstore;
}








/*
** helper functions for the higher level language to read the dotstore
*/
int GetDotX(DotStore *store, int index) { return store->GetDot(index)->x; }
int GetDotY(DotStore *store, int index) { return store->GetDot(index)->y; }
int GetDotLength(DotStore *store, int index) { return store->GetDot(index)->length; }
Dot *GetDot(DotStore *store, int index) { return store->GetDot(index); }
int GetNumDots(DotStore *store) { return store->GetNum(); }
void FreeDotStore(DotStore *store)  { delete store; }


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
	if(0)
	{
		const char *seq1="GATTACAATTAACTGATCGATCGTAGCTACATGCTGACTACTGACTGCATGCATGACTGCATGCATTGACTGACTGCATGACTGCATG";
		const char *seq2="AGCTCGATCGAGTCTCGAGTAG";
		int **pt=buildMappingTables(seq1,2);
	
		printf("%d - %d\n",(int)pt[0], (int)pt[1]);
	
		doComparison(pt, seq1, seq2, 2,10,2,6)->Dump();
	
		freeMappingTables(pt);
	}

	if(0)
	{
		const char *s1="GCGGGTACTGATATACTCATGATTATACCGCGCGGTTGTGTGAATTAATATCAACACCACAAAAGAGAGGAGGACTTCCTCTCTCTCTCTAACACCAATATATCCGGCCGGTTG";
		const char *s2="ATCGACGTATAGATTTTTCCACAGCGCCAAACTCTTCTATCACTCATGACTGACTGTGTCATGACTGATTATATATATCTCTCTTCTCATATATCATACT";
	
		printf("\nTEST1 should be 5\n");
		printf("matchAboveThreshold=%d\n",matchAboveThreshold(s1,12,s2,95,2,0,4) );
	
		printf("\nTEST2 should be 4\n");
		printf("matchAboveThreshold=%d\n",matchAboveThreshold(s1,24,s2,95,2,0,4) );
	
		char **results=convertSequence("GATACATTAAGCGC");
	
		printf(results[0]);
		printf("\n");
		printf(results[1]);
		printf("\n");
		printf(results[2]);
		printf("\n");
		
		printf("getTupleID('GGA')=%d\n",getTupleID("GGA",3));
	
		// k, thresh,wind
	// 	printf("matchAboveThreshold=%d\n",matchAboveThreshold(seq1+83,0,seq2+1,0,0,1,6) );
	// 	printf("%s\n%s\n",seq1+83,seq2+1);
	}
	
	const char *s1="GCGGGTACTGATATACTCATGATTATACCGCGCGGTTGTGTGAATTAATATCAACACCACAAAAGAGAGGAGGACTTCCTCTCTCTCTCTAACACCAATATATCCGGCCGGTTG";
	const char *s2="GCGGGTACTGATATACTCATGATTATACCGCGCGGTTGTGTGAATTAATATCAACACCACAAAAGAGAGGAGGACTTCCTCTCTCTCTCTAACACCAATATATCCGGCCGGTTG";
	//const char *s2="ATCGACGTATAGATTTTTCCACAGCGCCAAACTCTTCTATCACTCATGACTGACTGTGTCATGACTGATTATATATATCTCTCTTCTCATATATCATACT";
	
	DotStore *ds=makeDotComparisonByTranslation(s1, s2, /*int ktuplesize*/2, /*int window*/3, /*int mismatch*/0, /*int minmatch*/1);
	
	printf("%d matches\n",ds->GetNum());

	ds->Dump();

	printf("OLD\n");

	printf("%d matches\n",makeDotComparison(s1,s2,2,3,0,1)->GetNum());
	
}




}

