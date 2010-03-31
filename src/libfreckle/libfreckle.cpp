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

const char *Bases="ACGTN";							// N means UNKNOWN
const char *Aminos="ACDEFGHIKLMNPQRSTVWY-.";					// - = stop codon    . means UNKNOWN
const char *TranslateUniversal="KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV-Y-YSSSS-CWCLFLF";

/**
** \brief calculates the nth power of x.
** \details Uses tail recursion to calculate x to the power on n
** \param x the mantissa
** \param n the exponent
** \return x^n
*/
int ipow(int x, int n)
{
	return (!n)?1:(n&1)?x*ipow(x,n-1):ipow(x*x,n/2);
}

/**
** \brief calculates the nth power of x.
** \details Uses unsigned longs to calculate the nth power of x
** \param x the mantissa
** \param n the exponent
** \return x^n
*/
TupleStore uipow(TupleStore  base, TupleStore  exponent)
{
	if(exponent==0)
		return 1;
	exponent -= 1;
	return base * uipow(base, exponent);
}

/**
** \brief return the tuple id for len characters starting at tuple
** \details pass this a pointer to the base of a sequence string and a length.
** It will return the ktuple index of that sequence. It does this without building
** a lookup table to save memory.
** \param tuple The pointer into the string at the position where the tuple resides
** \param len How many characters are to be included in the tuple
** \param bases A pointer to the base pair character set. Default is "Bases" which indicates "ACGT" for Amino Acid sequences use "Aminos"
** \return returns an integer representing the tuple ID
*/
TupleID getTupleID(const char *tuple, int len, const char *bases)
{
	TupleID id=0;
	int index=0;

	for(int i=0; i<len; i++, tuple++)
	{
		index=strchr(bases, (int)*tuple)-bases;			// get the index in Bases of where this character appears
		if(!(index>=0 && index<(int)strlen(bases)))		// it should be within range
		{
			printf("index=%d bases=%d char=%c\n",index,(int)strlen(bases),*tuple);
			assert(0);
		}
		assert((index & BASE_MASK(bases)) == index);			// it shouldn't have any extra bits
		
		id=id * BASE_PAIRS(bases);		// 
	
		id += index;			// add the index bits onto the right
	}

	// assert that the tuple id is in a valid range
	assert(id>=0);
	assert(id<uipow(BASE_PAIRS(bases),len));

	// make sure it doesnt overflow
	assert(id+1);

	return id+1;				// 1 offset rather than 0
}

/**
** \brief Given a sequence string, build mapping tables "C" and "D"
** \details Given a sequence as a string, this function builds the mapping tables C and D as described in Huang and Zhangs paper
** \param sequence The entire sequence represented as a string of characters terminated by a NULL
** \param ktuplesize the size of the "ktuple word" in characters
** \param bases A pointer to the base pair character set. Use "Bases" or "Aminos" from the library
** \return returns a pointer to a table of two pointers. table[0] points to the location of C, and table[1] to D
*/
TupleStore **buildMappingTables( const char *sequence, int ktuplesize, const char *bases )
{
	int seqlen=strlen(sequence);
	assert(seqlen>0);

	TupleStore ktuplearraysize=uipow(BASE_PAIRS(bases),ktuplesize);
	int darraysize=seqlen-ktuplesize+1;
	
	// we allocate our arrays
	TupleStore *C=new TupleStore [ktuplearraysize];
	TupleStore *D=new TupleStore [darraysize];
	
	// lets zero our arrays
	memset(C, 0, sizeof(TupleStore)*ktuplearraysize);
	memset(D, 0, sizeof(TupleStore)*darraysize);

	// Initialise D
	const char *tuple=sequence;
	for(int i=0; i<darraysize; i++, tuple++)
		//read the ith tuple from the sequence. Put its tuple id into D
		D[i]=getTupleID(tuple,ktuplesize,bases);

	// Build the tables
	tuple=sequence;
	TupleStore cval=0;
	TupleStore dval=0;
	for(int i=0; i<darraysize; i++, tuple++)
	{
// 		printf("%d\n",i);
		// the index, i, of D is assigned to C[D[i]]. So if D[1]=3, then C[3]=1
		// of course all our C array offsets are 0 based, but the algorithm REQUIRES 1 based, because 0 is a terminator
		// if C[D[i]]!=0, then it already has an index. that value should be stored in D[i] before the C value is set
		dval=D[i];
// 		if(!dval<=ktuplearraysize)
//  			printf("%llu <= %llu\n",dval,ktuplearraysize);
		assert(dval<=ktuplearraysize);
		assert(dval>=1);
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
	
	TupleStore **pointertable=new TupleStore *[2];
	pointertable[0]=C;
	pointertable[1]=D;

	return pointertable;
}

/**
** \brief free the mapping tables as returned by buildMappingTables()
** \param tables the pointer two the two tables as returned by buildMappingTables()
** \return Nothing
** \see buildMappingTables
*/
void freeMappingTables(int **tables)
{
	delete tables[0];
	delete tables[1];
	delete tables;
}

/**
** \brief returns the sum of a buffer of ints
** \param buffer a pointer to the buffer
** \param length how many values to sum
** \return the sum
*/
int sum(int *buffer, int length)
{
	int sum=0;

	for(int i=0; i<length; i++)
		sum+=buffer[i];

	return sum;
}

/**
** \brief compute the matchlength of two subsequences
** \details given two sequences, seq1 and seq2, and positions in those sequences p1 and p2, compute how long the
** two sequences match for. this assumes that a length of 'k' matches already (the tuple size)
** \param seq1 The first sequence in which we are computing submatches
** \param p1 The position in the sequence where our initial match was found
** \param seq2 The second sequence
** \param p2 The position in the second sequence
** \param k The ktuple size used in the initial computation
** \param mismatch How many mismatches to allow during the scan per window size for it still to be considered "a match"
** \param window The size of the window for appraising mismatches.
*/
int matchAboveThreshold(const char *seq1, int p1, const char *seq2, int p2, int k, int mismatch, int window)
{
	assert(window>0);
	//assert(k>0);				//k=0 for methods 2 and 3 where we are using a coded ktuple to initiate these search locations
	int ringbuf[window];
	
	memset(ringbuf,0,sizeof(ringbuf));

	int matchlength=0;

	const char *s1=seq1+p1;
	const char *s2=seq2+p2;
	
	while((sum(ringbuf,window)<=mismatch) && *s1 && *s2 )
	{
		if(*s1=='.' || *s2=='.')						// unknowns are always a mismatch
		{
			ringbuf[matchlength++%window] = 1;
			s1++;
			s2++;
		}
		else 
			ringbuf[matchlength++%window] = *s1++==*s2++?0:1;		//compare s1 and s2 characters. assign 1 to the relevant 
	}	
	
	return matchlength-((sum(ringbuf,window)<=mismatch)?0:1);
		
}


/**
** \brief compare one sequence against the table constructed sequence
** \details once the tables "C" and "D" are created we can compare another (untabled) sequence against it using this function.
** This is the main workhorse function you will most probably use in libfreckle.
** \param tables the tables pointer returned from buildMappingTables()
** \param tablesequence the original sequence that was used to generate those tables
** \param newsequence the sequence to compare against the table sequence
** \param ktuplesize the size of the ktuple used during table construction
** \param window the window for appraising mismatches in matched subsequences
** \param mismatch how many characters per window can be allowed to mismatch for it still to be considered "matching"
** \param minmatch the minimum match length to store a dot for. This must be at least the size of the ktuple.
*/
DotStore *doComparison(unsigned int **tables, const char *tablesequence, const char *newsequence, int ktuplesize, int window, int mismatch, int minmatch, const char *bases )
{
	assert(mismatch<=window);
	assert(window>=ktuplesize);

	printf("doComparison(): tableseq=%d newseq=%d bases=%d\n",(int)strlen(tablesequence),(int)strlen(newsequence),(int)strlen(bases));

	unsigned int *C, *D;
	DotStore *dotstore=new DotStore();
	C=tables[0];
	D=tables[1];

	int newseqlen=(int)strlen(newsequence);
	assert(newseqlen>0);

	int darraysize=newseqlen-ktuplesize+1;

	// go through each k-tuple on the newsequence
	const char *tuple=newsequence;
	int tupleid=0;
	printf("%d\n",darraysize);
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
	unsigned int ***mappingtables=new unsigned int **[3];
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
				//int length=dot->length;

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

/*************************************************
** From here on is lbdot derivative code 
** but modified to fit our framework
 ************************************************/

#define DNAOrder       4

static char dnacode[256];
static char complment_base[256];

char *strrev( char *str)
{
int i, len=strlen(str);
char *p=new char [len+1];
    strcpy(p,str);
    for (i=0;i<len;i++) str[i]=p[len-i-1];
    str[len]=0; 
    delete p;
    return str;
};

void Init_code_tables()
{
	memset(dnacode,-1,sizeof (char) *256);
	dnacode['A']=dnacode['a']=0;
	dnacode['C']=dnacode['c']=1;
	dnacode['G']=dnacode['g']=2;
	dnacode['T']=dnacode['t']=3;

	memset(complment_base,'\0',sizeof (char) *256);
	complment_base['A']=complment_base['a']='T';
	complment_base['C']=complment_base['c']='G';
	complment_base['G']=complment_base['g']='C';
	complment_base['T']=complment_base['t']='A';
	complment_base['N']=complment_base['n']='N';

}

void EncodeNTSeq(const char *seq, int p1, int p2, int *c,int *d, int nm, int nMaxDNAKtup)
{          // c[i] contains last pos +1 of k_tuple No i
int i, j,m, L=p2-p1+1;
//char nt;
int v[64], sv;
	v[0]=1;
	for (i=1;i<nMaxDNAKtup+1;i++) v[i]=v[i-1]*DNAOrder;
	if(nm<1||nm>nMaxDNAKtup) nm=nMaxDNAKtup;
	sv=v[nm];
	for (j=0;j<=sv;j++) c[j]=0;

	d[0]=0;
	for(i=0;i<L;i++) {
		d[i+1]=dnacode[(int)seq[p1+i]];
	}
	L-=nm;

	for(i=1;i<=L;i++) {
		m=d[i];
		if(m>=0){
			for (j=1;j<nm;j++) {
				if (d[i+j]<0){
					m=-1; break;
				}
				m+=d[i+j]*v[j];
			}
		}
		d[i]=m;
	}

	for(i=1;i<=L;i++) {
		m=d[i];
		if(m<0) d[i]=0;
		else {
			d[i]=c[m]; c[m]=i;
		}
	}

}

int EncodeNTSeqConditional(const char *seq, int p1, int p2, int *c,int *d,int *cd, int nm, int maxHints, int nMaxDNAKtup)
{  //// c[i] contains last pos +1 of k_tuple No i
int i, j,m, L=p2-p1+1;
int v[64];
	v[0]=1;
	for (i=1;i<nMaxDNAKtup+1;i++) v[i]=v[i-1]*DNAOrder;
	if(nm<1||nm>nMaxDNAKtup) nm=nMaxDNAKtup;
int sv=v[nm];
int *stat=new int[sv+1];
	for (j=0;j<sv;j++) stat[j]=c[j]=0; 

	d[0]=0; cd[0]=0;
	for(i=0;i<L;i++) {
		d[i+1]=dnacode[(int)seq[p1+i]];
	}
	L-=nm;

	for(i=1;i<=L;i++) {
		m=d[i];
		if(m>=0){
			for (j=1;j<nm;j++) {
				if (d[i+j]<0){
					m=-1; break;
				}
				m+=d[i+j]*v[j];
			}
		}
		d[i]=m;
	}

	for(i=1;i<=L;i++) {
		m=d[i];
		if(m<0) d[i]=0;
		else {
			cd[i]=m;
			d[i]=c[m]; c[m]=i;
			(stat[m])++;
		}
	}

	int num=0;

	if(maxHints>100) {
		char txt[128];
		for(i=0;i<=sv;i++) {
			if(stat[i]>maxHints) {
				strncpy(txt,seq+c[i]-1,nm);txt[nm]=0;
				printf("%s,code %d repeats %d times\n",txt,i, stat[i]);
				c[i]=0;
				num++;
			}
		}

	}

	delete stat;

	return num;
}

int GetNtCode(const char *seq, int ktup, int intval, const int *v)
{
	int i, m=0, c;
	for (i=0;i<ktup;i++) {
		c=dnacode[(int)seq[i*intval]]*v[i];
		if (c<0) return -1;
		m+=c;
	}
	return m;
};

void ComplementSeq(char  *a)
{
	while (*a) { 
	  *a=complment_base[(int)*a];
	  a++;
	}
}

char *RCseq(char *a)
{
	strrev(a);
	ComplementSeq(a);
	return a;

}

// This is the obfuscated lbdot original
DotStore **DoFastComparison(char *Seq1, char *Seq2, int SeqLen1, int SeqLen2, int CompWind,int CompMism, int nMaxRepeatKtup, int nMaxDNAKtup)
{
// printf("initcodetables, %d, %d, %d, %d (%d,%d)\n",CompWind, CompMism, nMaxRepeatKtup, nMaxDNAKtup,SeqLen1,SeqLen2);
// 
// printf("start seq 1: %s\n",strndup(Seq1,80));
// printf("start seq 2: %s\n",strndup(Seq2,80));
// printf("end seq 1: %s\n",strndup(Seq1+strlen(Seq1)-80,80));
// printf("end seq 2: %s\n",strndup(Seq2+strlen(Seq2)-80,80));
// 
// printf("extra characters seq1:\n");
// for(char *c=Seq1; *c; c++)
// 	if(*c != 'A' && *c != 'C' && *c != 'G' && *c != 'T')
// 		printf("%c",*c);
// printf("\n\n");
// printf("extra characters seq2:\n");
// for(char *c=Seq2; *c; c++)
// 	if(*c != 'A' && *c != 'C' && *c != 'G' && *c != 'T')
// 		printf("%c",*c);
// printf("\n\n");

DotStore **result=new DotStore *[2];

Init_code_tables();

DotStore *PlusDotArray=new DotStore();
DotStore *MinusDotArray=new DotStore();

int x,i,ix,j,cmpNum, ic, ct, ctt, nBreak;
int CompUnit=CompWind;
int CompErr=CompMism;
int CompKtup;
int Length1=SeqLen1;
int Length2=SeqLen2;
bool bRCSeq=false, ok=true;

int Reject1=0, Reject2=0;
int nQualified=0;
Dot cp;
char *s1= Seq1;
char *s2= Seq2;
int ONE=1, *pc;
int *d1, *c1;
int v[64], pos, scs, dp=-CompErr;
char *sc=new char [CompUnit+2];

	if(Seq1==Seq2){
		s2=new char [Length2+2];
		strcpy(s2, Seq2);
	}

	cmpNum=2;

	v[0]=1;
	CompKtup=(CompUnit<nMaxDNAKtup)?CompUnit:nMaxDNAKtup;
	for (i=1;i<=nMaxDNAKtup;i++) v[i]=v[i-1]*4;

// 	printf("Dotplot Fast Method, win=%d, mismatch=%d, ktup=%d\n",
// 			CompWind,CompMism,CompKtup);

	int pm=v[CompKtup];
	c1=new int [pm+2];
	d1=new int [Length1+2];
int *cd=new int [Length1+2];
	for (i=0;i<=Length1; i++) cd[i]=0;

	if(nMaxRepeatKtup>200){// may significantly decrease computing for repeatitive seq
		(void)EncodeNTSeqConditional(s1,0,Length1-1, c1,d1,cd,CompKtup, nMaxRepeatKtup,nMaxDNAKtup);
	} else {
		EncodeNTSeq(s1,0,Length1-1, c1,d1, CompKtup, nMaxDNAKtup);
		/// make sure c1[cd[]]>0
		j=GetNtCode(s1,CompKtup,1, v);
		for (i=0;i<Length1;i++) cd[i]=j;
	}

// 	printf("Comparing...\n");
	for (ic=0;ic<cmpNum&&ok;ic++)  {
		Reject1=Reject2=0;
		nQualified=0;
		if(ic==1) {
			bRCSeq=true;
			RCseq(s2);
		}

		pc=(Seq1==Seq2&&!bRCSeq)?&j:&ONE;


		int dd=Length2-CompKtup;
		int ddt=dd/100+1;

		if(Length2<pm*2){//bUseLessMem
		////////when Length2>pm*2 the next method may run faster 
		////////due to repetitive instructions with lookup table
//////////////////////////////////////////////////////////////////////
/////////////Here is the different block from above method/////+ j=d2[j] }
			for(j=1;j<dd&&ok;j++){
				i=GetNtCode(s2+j-1,CompKtup,1, v);
				if (i<0) continue;
////////////////////////////////////////////////////////////
				ix=c1[i];
				while (ix>=*pc) {
					if((j>1&&ix>1)&&s2[j-2]==s1[ix-2]){
						// if previous bases match, ignore current dot 
						//// however, if previous pair was ignored due to high repeats,
						// don't give up the current dot
						if(c1[cd[ix-1]]>0||
							(pc==&j&&j==ix)){/// same seq on diagnal
							ix=d1[ix];
							Reject1++;
							continue;
						}
					}
					ctt=Length2-j+1;
					ct=Length1-ix+1;
					if(ctt>ct) ctt=ct;
					ct=CompKtup;
					while(ct<ctt&&(s2[ct+j-1]==s1[ct+ix-1])) 	ct++;


					if(CompErr>0){
						scs=0;
						memset(sc,0,CompUnit); //// match=0 mimatch=-1
						pos=ct%CompUnit;
						while(ct<ctt&&scs>=dp){
							scs-=sc[pos];
							if(s2[ct+j-1]==s1[ct+ix-1]) sc[pos]=0;
							else sc[pos]=-1;
							scs+=sc[pos];

//							pos=(pos+1)%CompUnit;/// this is slow
							pos++; 
							if(pos>=CompUnit) pos=0;
							ct++;
						}
					}
					
					if(ct<CompUnit){
						ix=d1[ix];
						Reject2++;
						continue; /// not long enough
					}
					nQualified++;

					if(CompErr>0){/// try to extend ct as far as possible
						nBreak=0;
						while(ct<ctt&&nBreak<CompUnit){
							scs-=sc[pos];
							if(s2[ct+j-1]==s1[ct+ix-1]) sc[pos]=0;
							else sc[pos]=-1;
							scs+=sc[pos];

							pos++;
							if(pos>=CompUnit) pos=0;
							ct++;

							if(scs<dp) nBreak++;
							else nBreak=0;
						}
						ct-=nBreak;

					}
					
					cp.x=ix-1;cp.y=j-1;cp.length=ct;

					if(bRCSeq) {
						MinusDotArray->AddDot(cp.x,cp.y,cp.length);
					} else {
						PlusDotArray->AddDot(cp.x,cp.y,cp.length);
							if(*pc>1){/// Add mirror point
								x=cp.x; cp.x=cp.y; cp.y=x;
								PlusDotArray->AddDot(cp.x,cp.y,cp.length);
								x=cp.x; cp.x=cp.y; cp.y=x;
							}
					}

					ix=d1[ix];
				}
//				if (j%ddt==0){
//					printf("%d%%", ic*50+MulDiv(j,100,dd*cmpNum));
//				}
			}
		} else {
			int *c2=NULL;
			int *d2=NULL; 

			if(Seq1==Seq2&&!bRCSeq){
				c2=c1; d2=d1;
				pc=&j;
			} else {
				pc=&ONE;
				c2=new int [pm+2];d2=new int [Length2+2];
				EncodeNTSeq(s2,0,Length2-1, c2,d2, CompKtup, nMaxDNAKtup);
			}


			ddt=pm/100+1;
//////////////////////////////////////////////////////////////////////
/////////////Here is the different block from above method/////+ j=d2[j] }
			for (i=0;i<pm&&ok;i++){
				if(c1[i]<1) continue; /// ignore if the other seq no such k-tuple
				j=c2[i];
				while(j>0){
////////////////////////////////////////////////////////////
					ix=c1[i];
					while (ix>=*pc) {
						if((j>1&&ix>1)&&s2[j-2]==s1[ix-2]) {
						// if previous bases match, ignore current dot 
						//// however, if previous pair was ignored due to high repeats,
						// don't give up the current dot
							if(c1[cd[ix-1]]>0||
								(pc==&j&&j==ix)){/// same seq on diagnal
								ix=d1[ix];
								Reject1++;
								continue;
							}
						}
						ctt=Length2-j+1;
						ct=Length1-ix+1;
						if(ctt>ct) ctt=ct;
						ct=CompKtup;
						while(ct<ctt&&(s2[ct+j-1]==s1[ct+ix-1])) 	ct++;

						if(CompErr>0){
							scs=0;
							memset(sc,0,CompUnit); //// match=0 mimatch=-1
							pos=ct%CompUnit;
							while(ct<ctt&&scs>=dp){
								scs-=sc[pos];
								if(s2[ct+j-1]==s1[ct+ix-1]) sc[pos]=0;
								else sc[pos]=-1;
								scs+=sc[pos];

//							pos=(pos+1)%CompUnit;/// this is slow
								pos++; 
								if(pos>=CompUnit) pos=0;
								ct++;
							}
						}
					
						if(ct<CompUnit){
							ix=d1[ix];
							Reject2++;
							continue; /// not long enough
						}
						nQualified++;

						if(CompErr>0){/// try to extend ct as far as possible
							nBreak=0;
							while(ct<ctt&&nBreak<CompUnit){
								scs-=sc[pos];
								if(s2[ct+j-1]==s1[ct+ix-1]) sc[pos]=0;
								else sc[pos]=-1;
								scs+=sc[pos];

								pos++;
								if(pos>=CompUnit) pos=0;
								ct++;

								if(scs<dp) nBreak++;
								else nBreak=0;
							}
							ct-=nBreak;
						}
					
						cp.x=ix-1;cp.y=j-1;cp.length=ct;

						if(bRCSeq) {
							MinusDotArray->AddDot(cp.x,cp.y,cp.length);
						} else {
							PlusDotArray->AddDot(cp.x,cp.y,cp.length);
							if(*pc>1){/// Add mirror point
								x=cp.x; cp.x=cp.y; cp.y=x;
								PlusDotArray->AddDot(cp.x,cp.y,cp.length);
								x=cp.x; cp.x=cp.y; cp.y=x;
							}
						}

						ix=d1[ix];
					}
//					if (i%ddt==0){
//						printf("%d%%", ic*50+MulDiv(i,100,pm*cmpNum));
//					}
					j=d2[j];
				}
			}
			if(c2!=c1) delete c2;
			if(d2!=d1) delete d2;
		}
		if(bRCSeq) 	RCseq(s2);
	}
	delete c1; delete d1;delete sc; delete cd;

	if(Seq1==Seq2)delete s2;

// 	delete MinusDotArray;
	
	//PlusDotArray->Interpolate(CompWind);
	//MinusDotArray->Interpolate(CompWind);

// 	return PlusDotArray;
	result[0]=PlusDotArray;
	result[1]=MinusDotArray;
	return result;
/*
	return CompKtup;*/
};

/*
** helper functions for the higher level language to read the dotstore
*/
DotStore *NewDotStore() { return new DotStore(); }
void DelDotStore(DotStore *store) { delete store; } 
void DotStoreAddDot(DotStore *store, int x, int y, int len) { store->AddDot(x,y,len); }
int DotStoreGetDotX(DotStore *store, int index) { return store->GetDot(index)->x; }
int DotStoreGetDotY(DotStore *store, int index) { return store->GetDot(index)->y; }
int DotStoreGetDotLength(DotStore *store, int index) { return store->GetDot(index)->length; }
Dot *DotStoreGetDot(DotStore *store, int index) { return store->GetDot(index); }
int DotStoreGetNumDots(DotStore *store) { return store->GetNum(); }
void DotStoreCreateIndex(DotStore *store) { store->CreateIndex(); }
void DotStoreDestroyIndex(DotStore *store) { store->DestroyIndex(); }
int *DotStoreToBuffer(DotStore *store) { return store->ToBuffer(); }
void DotStoreFromBuffer(DotStore *store, int *buffer) { store->FromBuffer(buffer); }
int DotStoreBufferSize(DotStore *store, int *buffer) { return store->BufferSize(buffer); }
void FreeIntBuffer(int *buffer) { assert(buffer); delete buffer; }
DotStore *DotStoreFilter(DotStore *store, int minlen) { return store->Filter(minlen); } 
void DotStoreInterpolate(DotStore *store, int window) { store->Interpolate(window); }

void DotStoreSetMaxX(DotStore *store, int max) {store->SetMaxX(max);}
void DotStoreSetMaxY(DotStore *store, int max) {store->SetMaxY(max);}
int DotStoreGetMaxX(DotStore *store) {return store->GetMaxX();}
int DotStoreGetMaxY(DotStore *store) {return store->GetMaxY();}

// conservation helper functions
Dot *DotStoreGetIndexLongestMatchingRowDot(DotStore *store, int x) { return store->GetIndexLongestMatchingRowDot(x); }
Dot *DotStoreGetIndexLongestMatchingColumnDot(DotStore *store, int y) { return store->GetIndexLongestMatchingColumnDot(y); }

/*
** helper function to interface with the dotgrid
*/
DotGrid *NewDotGrid() { return new DotGrid(); }
void DelDotGrid(DotGrid *grid) { delete grid; }
void DotGridCreate(DotGrid *grid, int width, int height) { grid->Create(width, height); }
void DotGridSetPoint(DotGrid *grid, int x, int y, int point) { grid->SetPoint(x,y,point); }
int DotGridGetPoint(DotGrid *grid, int x, int y) { return grid->GetPoint(x,y); }
int DotGridGetData(DotGrid *grid, int pos) { return grid->GetData(pos); }
int DotGridWidth(DotGrid *grid) { return grid->GetWidth(); }
int DotGridHeight(DotGrid *grid) { return grid->GetHeight(); }
int DotGridGetSize(DotGrid *grid) { return grid->GetSize(); }
int DotGridGetMax(DotGrid *grid) { return grid->GetMax(); }
int DotGridGetMin(DotGrid *grid) { return grid->GetMin(); }
unsigned char *DotGridToString(DotGrid *grid) { return grid->ToString(); }
void FreeString(unsigned char *string) { assert(string); delete string; }
void DotGridCalculate(DotGrid *grid, DotStore *source, double x1, double y1, double x2, double y2, double scale, int window)
{
	grid->CalculateGrid(source, x1, y1, x2, y2, scale, window);
}
void DotGridAddInplace(DotGrid *source, DotGrid *add) { source->AddInplace(add); }
void DotGridFlipInplace(DotGrid *grid) { grid->FlipInplace(); }
int *DotGridCalculateHistogram(DotGrid *grid) {return grid->CalculateHistogram();}

// unsigned char *DotStoreImageToString(DotStore *store, int xseqsize, int yseqsize, int longest, int window)
// {
// 	store->CreateIndex();
// 	int *storage=store->CalculateAverageGrid(xseqsize,yseqsize,longest,window);
// 	unsigned char *result=store->GridToString();
// 	return result;
// }

void DumpDotStore(DotStore *store) { store->Dump(); }

#ifdef DEBUG
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
	makeDotComparison(s1,s2,2,3,0,1)->Dump();
// 	printf("%d matches\n",makeDotComparison(s1,s2,2,3,0,1)->GetNum());
	
}
#endif



}

