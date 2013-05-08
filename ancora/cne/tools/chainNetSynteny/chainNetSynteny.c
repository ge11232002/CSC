/* chainNetSynteny - Make synteny blocks from alignment nets. 
 * This program was originally written by Par Engstrom, BCCS, University of Bergen, Norway in 2006
 * for a paper published in 2007 (http://www.genome.org/cgi/content/abstract/17/12/1898).
 * It heavily uses the UCSC Genome Browser library (http://hgdownload.cse.ucsc.edu/admin/jksrc.zip).
 * The synteny block detection pipeline, of which this program is a component, is described in the paper.
 */

#include "common.h"
#include "linefile.h"
#include "localmem.h"
#include "hash.h"
#include "options.h"
#include "dnautil.h"
#include "rbTree.h"
#include "chain.h"
#include "chainNet.h"
#include "portable.h"
#include "bed.h"


struct optionSpec optionSpecs[] =
{
    {"frag", OPTION_BOOLEAN },
    {"minAli", OPTION_INT },
    {"breakingAli", OPTION_INT },
    {"minGapShown", OPTION_INT },
    {"tTrackOpt", OPTION_STRING },
    {"qTrackOpt", OPTION_STRING },
    {"longName", OPTION_BOOLEAN },
    {NULL, 0}
};

boolean queryIsFragmented = 0; /* Whether or not to allow superblocks to across query sequences */
boolean writeLongName = 0; /* Whether or not to write long names in bed file */
int minAli = 5000;	/* Minimum gap size to fill. */
int minBreakingAli = 5000; /* Minimum fill to record. */
int minGapShown = 0; /* Minimum gap size to show in BED output (0=hide all gaps) */
char *tTrackOptions = NULL;
char *qTrackOptions = NULL;

void usage()
/* Explain usage and exit. */
{
    errAbort(
	"chainNetSynteny - Make synteny blocks from alignment nets\n"
	"usage:\n"
	"   chainNetSynteny in.chain target.net query.net target.bed query.bed\n"
	"where:\n"
	"   in.chain is the stitched chain file for reciprocal nets\n"
	"   target.net is the one-way nets for the target genome\n"
	"   query.net is the one-way nets for query genome\n"
        "   target.bed is the output for the target genome\n"
        "   query.bed is the output for the query genome\n"
	"options:\n"
	"   -frag             - Query sequences are fragmented.\n"
        "                       Create synteny blocks across query sequences.\n"
	"   -minAli=N         - Default %d.\n"
	"   -breakingAli=N    - Default %d.\n"
	"   -minGapShown=N    - Min gap size to show in bed output.\n"
        "                       0=hide all gaps. Default %d.\n"
	"   -tTrackOpt=string - Track line options for target.bed.\n"
	"   -qTrackOpt=string - Track line options for query.bed.\n"
	"   -longName         - Write long names in bed file.\n",
	minAli, minBreakingAli, minGapShown);
}


struct alignedRegion
{
    int start, end;
};


struct syntenyBlock
{
    struct syntenyBlock *next;
    struct chain *chainList;
};


char *makeAbbreviatedLocString(char *chrom, int pos)
/* Make a string of format chr1:54200k */
{
    int l = strlen(chrom);
    char *s = needMem(l+9); /* space for chrom + ':' + 6 digits + 'k' */
    pos = pos / 1000;
    assert(pos < 1000000);  /* make sure pos is <= 6 digits long (i.e. < 1 gigabase) */
    sprintf(s, "%s:%dk", chrom, pos);
    return s;
}


char *makeLongTargetLocString(struct chain *chain)
/* Make a string of format chr1:2000-2100,chr2:300-500 for target coordinates of chain list */
{
    char *s1 = NULL;
    char *s2;

    for(; chain; chain = chain->next)
    {
	assert(chain->tStart < 999999999); /* make sure pos is <= 9 digits long (i.e. < 1 gigabase) */
	assert(chain->tEnd < 999999999);   /* make sure pos is <= 9 digits long (i.e. < 1 gigabase) */
	if(s1) {
	    s2 = needMem(strlen(s1) + 1 + strlen(chain->tName) + 1 + 9 + 1 + 9 + 1); 
	    sprintf(s2, "%s,%s:%d-%d", s1, chain->tName, chain->tStart+1, chain->tEnd);
	    freeMem(s1);
	    s1 = s2;
	}
	else {
	    s1 = needMem(strlen(chain->tName) + 1 + 9 + 1 + 9 + 1); 
	    sprintf(s1, "%s:%d-%d", chain->tName, chain->tStart+1, chain->tEnd);
	}
    }
   
    return s1;
}


char *makeLongQueryLocString(struct chain *chain)
/* Make a string of format chr1:2000-2100,chr2:300-500 for query coordinates of chain list */
{
    char *s1 = NULL;
    char *s2;
    int start, end;

    for(; chain; chain = chain->next)
    {
	if(chain->qStrand == '+') {
	    start = chain->qStart;
	    end = chain->qEnd;
	}
	else {
	    start = chain->qSize - chain->qEnd;
	    end = chain->qSize - chain->qStart;
	}
	assert(start < 999999999); /* make sure pos is <= 9 digits long (i.e. < 1 gigabase) */
	assert(end < 999999999);   /* make sure pos is <= 9 digits long (i.e. < 1 gigabase) */
	if(s1) {
	    s2 = needMem(strlen(s1) + 1 + strlen(chain->qName) + 1 + 9 + 1 + 9 + 1); 
	    sprintf(s2, "%s,%s:%d-%d", s1, chain->qName, start+1, end);
	    freeMem(s1);
	    s1 = s2;
	}
	else {
	    s1 = needMem(strlen(chain->qName) + 1 + 9 + 1 + 9 + 1); 
	    sprintf(s1, "%s:%d-%d", chain->qName, start+1, end);
	}
    }
   
    return s1;
}


int chainCmpQueryPlusStrand(const void *va, const void *vb)
/* Compare to sort based on query chrom and position, accounting for strand. */
{
const struct chain *a = *((struct chain **)va);
const struct chain *b = *((struct chain **)vb);
int dif;                                                                        

dif = strcmp(a->qName, b->qName);                                               
if (dif == 0)                                                                   
    dif =
	(a->qStrand == '+' ? a->qStart : a->qSize - a->qEnd) -
	(b->qStrand == '+' ? b->qStart : b->qSize - b->qEnd);
return dif;                       
}


int alignedRegionCmp(void *va, void *vb)
/* Return -1 if a before b,  0 if a and b overlap,
 * and 1 if a after b. */
{
    struct alignedRegion *a = va;
    struct alignedRegion *b = vb;
    if (a->end <= b->start)
	return -1;
    else if (b->end <= a->start)
	return 1;
    else
	return 0;
}


void rFindAlignedRegions(struct cnFill *fillList, struct rbTree *ari)
/* Recursively find aligned regions in net and add them to the supplied index */
{
    struct cnFill *fill, *gap;
    struct alignedRegion *ar;
    int tStart;

    /* Iterate over all fills on this level */
    for (fill = fillList; fill != NULL; fill = fill->next)
    {
	tStart = fill->tStart;

	/* Iterate over all gaps in the fill to process child fills and
	   compute which regions of the fill are outside gaps */
	for (gap = fill->children; gap != NULL; gap = gap->next)
	{
	    /* Recurse to process fills inside the gap, if any */
	    if(gap->children)
		rFindAlignedRegions(gap->children, ari);

	    /* Add an aligned region for the current fill
              (from end of previous gap up to start of current gap) */
	    AllocVar(ar);
	    ar->start = tStart;
	    ar->end = gap->tStart;
	    rbTreeAdd(ari, ar);
	    //printf(" %d-%d %s\n", ar->start, ar->end, fill->qName);

	    /* Move to end of gap */
	    tStart = gap->tStart + gap->tSize;
	}

	/* Add the final aligned region for the current fill */
	AllocVar(ar);
	ar->start = tStart;
	ar->end = fill->tStart+fill->tSize;
	rbTreeAdd(ari, ar);
	//printf(" %d-%d %s\n", ar->start, ar->end, fill->qName);
    }
}


struct hash *buildAlignedRegionIndex(char *netFile)
{
    struct lineFile *lf = lineFileOpen(netFile, TRUE);
    struct hash *genomeAri = newHash(0);
    struct rbTree *chromAri;
    struct chainNet *net;

    while ((net = chainNetRead(lf)) != NULL)
    {
	verbose(2, " finding aligned regions for %s...\n", net->name);
	
	/* Create index for chromosome */
	chromAri = rbTreeNew(alignedRegionCmp);
	rFindAlignedRegions(net->fillList, chromAri);

	/* Add chrom index to genome index */
	hashAddUnique(genomeAri, net->name, chromAri);

	/* Free net structure */
	chainNetFree(&net);
    }

    return genomeAri;
}

    
int countAlignedBasesInRegion(struct rbTree *ari, int start, int end)
{
    struct slRef *arListItem, *arList;
    struct alignedRegion *ar, startRegion, endRegion;
    int count = 0;

    startRegion.start = start;
    startRegion.end = start+1;
    endRegion.start = end-1;
    endRegion.end = end;
    arList = rbTreeItemsInRange(ari, &startRegion, &endRegion);

    for(arListItem = arList; arListItem != NULL; arListItem = arListItem->next)
    {
	ar = arListItem->val;
	count += (min(end, ar->end) -  max(start, ar->start));
	/*if(count > minFill) 
	{
	    slFree(&arList);
	    return 0;
	    }*/
	    
    }
    
    slFreeList(&arList);
    return count;
}


int countAlignedBasesInChain(struct chain *chain)
{
    struct cBlock *cBlock;
    int ali = 0;
    for(cBlock = chain->blockList; cBlock != NULL; cBlock = cBlock->next) {
	ali += (cBlock->tEnd - cBlock->tStart);
    }
    return ali;
}


void breakChainAtFilledGaps(struct chain **pChain, struct rbTree *tAri, struct rbTree *qAri, struct chain **chainList)
/* Find synteny blocks in a chain and free the original chain */
{
    struct chain *chain = *pChain;
    struct chain *subChain, *chainToFree;
    struct cBlock *cBlock1, *cBlock2, *syntenyStartBlock;
    int d1, d2, ali;

    syntenyStartBlock = cBlock1 = chain->blockList;
    ali = cBlock1->tEnd - cBlock1->tStart;
    for(cBlock2 = cBlock1->next; cBlock2 != NULL; cBlock2 = cBlock2->next) 
    {
	d1 = cBlock2->tStart - cBlock1->tEnd;
	d2 = cBlock2->qStart - cBlock1->qEnd;
	if((d1 >= minBreakingAli &&
	    countAlignedBasesInRegion(tAri, cBlock1->tEnd, cBlock2->tStart) >= minBreakingAli)
	   ||
	   (d2 >= minBreakingAli &&
	    (chain->qStrand == '+' ?
	     countAlignedBasesInRegion(qAri, cBlock1->qEnd, cBlock2->qStart) :
	     countAlignedBasesInRegion(qAri, chain->qSize - cBlock2->qStart, chain->qSize - cBlock1->qEnd))
	    >= minBreakingAli))
	{
	    if(ali >= minAli)
	    {
		chainSubsetOnT(chain, syntenyStartBlock->tStart, cBlock1->tEnd, &subChain, &chainToFree);
		slAddHead(chainList, subChain);
		/*if(syntenyStartBlock != chain->blockList)
		    fprintf(breakOut, "%s\t%d\t%d\n", chain->tName, syntenyStartBlock->tStart, syntenyStartBlock->tStart+1);
		fprintf(breakOut, "%s\t%d\t%d\n", chain->tName, cBlock1->tEnd-1, cBlock1->tEnd);*/
	    }
	    syntenyStartBlock = cBlock2;
	    ali = 0;
	}
	cBlock1 = cBlock2;
	ali += cBlock1->tEnd - cBlock1->tStart;
    }

    if(ali >= minAli) {
	chainSubsetOnT(chain, syntenyStartBlock->tStart, cBlock1->tEnd, &subChain, &chainToFree);
	slAddHead(chainList, subChain);
	/*if(syntenyStartBlock != chain->blockList)
	    fprintf(breakOut, "%s\t%d\t%d\n", chain->tName, syntenyStartBlock->tStart, syntenyStartBlock->tStart+1);*/
    }

    /* The only case where the original chain should not be freed is if it is identical to the
       chain was made into a block */
    if(chainToFree != NULL)
	chainFree(pChain);
}


void removeChainOverlaps(struct chain **chainList)
{
    struct chain **prevNextPointer;
    struct chain *chain1, *chain2, *subChain, *chainToFree;
    char *chrom;
    int minStart, maxEnd;
    int chain1PlusStart, chain2PlusStart, chain1PlusEnd, chain2PlusEnd;   

    /* Sort list on target position and traverse list to find overlaps in target */
    slSort(chainList, chainCmpTarget);
    minStart = -1;
    chrom = "";
    prevNextPointer = chainList;
    while(*prevNextPointer)
    {
	/* Set chain1 to the current position in the list */
	chain1 = *prevNextPointer;
        /* Set chain2 to the following chain on the same chrom, or NULL */
	chain2 = (chain1->next && sameString(chain1->tName, chain1->next->tName)) ?  chain1->next : NULL;
	/* Reset minStart if the previous chain was on  a different chrom */
	if(!sameString(chrom, chain1->tName)) {
	    chrom = cloneString(chain1->tName);
	    minStart = -1;
	}
	/* Set minStart and maxEnd to bounds of (truncated) chain */
	minStart = max(minStart, chain1->tStart);
	maxEnd = chain2==NULL ? chain1->tEnd : min(chain1->tEnd, chain2->tStart);
	/* Shall we truncate the chain? */
	if(chain1->tStart < minStart || chain1->tEnd > maxEnd)
	{ /* Yes: truncate */
	    verbose(2, "Truncating chain %d from %s:%d-%d to %d-%d\n",
		    chain1->id, chain1->tName, chain1->tStart, chain1->tEnd, minStart, maxEnd);
	    /* Generate truncated chain and insert it into the list if it is large enough */
	    if(maxEnd - minStart >= minAli) /* Only attempt truncation if the result may be large enough */
	    {
		chainSubsetOnT(chain1, minStart, maxEnd, &subChain, &chainToFree);
		assert(chainToFree != NULL); /* Make sure we got a truncated copy */
		if(countAlignedBasesInChain(subChain) >= minAli)
		{ /* Insert before chain1 */
		    subChain->next = chain1;
		    *prevNextPointer = subChain;
		    prevNextPointer = &(subChain->next);
		    verbose(2, " (inserted subchain)\n");
		}
	    }
	    /* Does chain1 extend past chain2? */
	    if(chain2 != NULL && chain1->tEnd > chain2->tEnd)
	    { /* Yes: unlink and free chain2 */
		verbose(2, "Removing chain %d at %s:%d-%d\n",
			chain2->id, chain2->tName, chain2->tStart, chain2->tEnd);
		minStart = chain2->tEnd;
		chain1->next = chain2->next;
		chainFree(&chain2);
	    }
	    else
	    { /* No: unlink and free chain1 */
		minStart = chain1->tEnd;
		*prevNextPointer = chain1->next;
		chainFree(&chain1);
	    }
	}
	else
	{ /* No: don't trucate, just move forward by updating prevNextPointer */
	    prevNextPointer = &(chain1->next);
	}
    }

    /* Repeat for query genome, i.e.
     * sort list on query position and traverse to find overlaps in query */
    slSort(chainList, chainCmpQueryPlusStrand); 
    minStart = -1;
    chrom = "";
    prevNextPointer = chainList;
    while(*prevNextPointer)
    {
	/* Set chain1 to the current position in the list */
	chain1 = *prevNextPointer;
	if(chain1->qStrand == '+')
	{
	    chain1PlusStart = chain1->qStart;
	    chain1PlusEnd = chain1->qEnd;
	}
	else
	{
	    chain1PlusStart = chain1->qSize - chain1->qEnd;
	    chain1PlusEnd = chain1->qSize - chain1->qStart;
	}
        /* Set chain2 to the following chain on the same chrom, or NULL */
	if(chain1->next && sameString(chain1->qName, chain1->next->qName))
	{
	    chain2 = chain1->next;
	    if(chain2->qStrand == '+')
	    {
		chain2PlusStart = chain2->qStart;
		chain2PlusEnd = chain2->qEnd;
	    }
	    else
	    {
		chain2PlusStart = chain2->qSize - chain2->qEnd;
		chain2PlusEnd = chain2->qSize - chain2->qStart;
	    }
	}
	else
	{
	    chain2 = NULL;
	    chain2PlusStart = chain2PlusEnd = 0; /* To avoid compiler warnings */
	}
	/* Reset minStart if the previous chain was on a different chrom */
	if(!sameString(chrom, chain1->qName)) {
	    chrom = cloneString(chain1->qName);
	    minStart = -1;
	}
	/* Set minStart and maxEnd to bounds of (truncated) chain */
	minStart = max(minStart, chain1PlusStart);
	maxEnd = chain2==NULL ? chain1PlusEnd : min(chain1PlusEnd, chain2PlusStart);
	/* Shall we truncate the chain? */
	if(chain1PlusStart < minStart || chain1PlusEnd > maxEnd)
	{ /* Yes: truncate */
	    verbose(2, "Truncating chain %d from %s:%d-%d to %d-%d\n",
		    chain1->id, chain1->qName, chain1PlusStart, chain1PlusEnd, minStart, maxEnd);
	    /* Generate truncated chain and insert it into the list if it is large enough */
	    if(maxEnd - minStart >= minAli) /* Only attempt truncation if the result may be large enough */
	    {
		if(chain1->qStrand == '+')
		    chainSubsetOnQ(chain1, minStart, maxEnd, &subChain, &chainToFree);
		else
		    chainSubsetOnQ(chain1, chain1->qSize - maxEnd, chain1->qSize - minStart, &subChain, &chainToFree);
		assert(chainToFree != NULL); /* Make sure we got a truncated copy */
		if(countAlignedBasesInChain(subChain) >= minAli)
		{ /* Insert before chain1 */
		    subChain->next = chain1;
		    *prevNextPointer = subChain;
		    prevNextPointer = &(subChain->next);
		    verbose(2, " (inserted subchain)\n");
		}
	    }
	    /* Does chain1 extend past chain2? */
	    if(chain2 != NULL && chain1PlusEnd > chain2PlusEnd)
	    { /* Yes: unlink and free chain2 */
		verbose(2, "Removing chain %d at %s:%d-%d\n",
			chain2->id, chain2->qName, chain2PlusStart, chain2PlusEnd);
		minStart = chain2PlusEnd;
		chain1->next = chain2->next;
		chainFree(&chain2);
	    }
	    else
	    { /* No: unlink and free chain1 */
		minStart = chain1PlusEnd;
		*prevNextPointer = chain1->next;
		chainFree(&chain1);
	    }
	}
	else
	{ /* No: don't trucate, just move forward by updating prevNextPointer */
	    prevNextPointer = &(chain1->next);
	}
    }
}


struct syntenyBlock *createSyntenyBlocks(struct chain **chainList, struct hash *tGenomeAri, struct hash *qGenomeAri)
{
    struct rbTree *chromAri;
    struct syntenyBlock *sb, *sbList = NULL;
    struct chain *chain, *next;
    int d1, d2, ali;
    int chainCount = 0;

    slSort(chainList, chainCmpTarget); /* Sort the list of chains on target position */
    chain = *chainList; /* Get the first chain in the list */
    *chainList = NULL;  /* NULL the list pointer because the list will be emptied */

    /* Loop over all chains in the list and group them into synteny blocks */
    while(chain)
    {
	/* Create a new synteny block and assign it the entire remaining chain list */
	AllocVar(sb);
	sb->chainList = chain;
	/* Add the synteny block to the list of synteny blocks */
	slAddHead(&sbList, sb);
        /* Loop through the chain list until we find a chain that is not colinear with the next */
	for(;;)
	{
	    chainCount++;
	    next = chain->next;
	    /* Break if there are no more chains or next chain is on a different target chrom */
	    if(next == NULL || !sameString(chain->tName, next->tName))
		break;
	    /* Is next chain on the same query chrom? */
	    if(sameString(chain->qName, next->qName))
	    { /* Yes */
		/* Break if different query strand or non-increasing order on query sequence */
		if(chain->qStrand != next->qStrand || chain->qEnd > next->qStart)
		    break;
		/* Break if number of aligned target bases between chains is >= minBreakingAli */
		if(next->tStart - chain->tEnd >= minBreakingAli)
		{
		    chromAri = hashMustFindVal(tGenomeAri, chain->tName);
		    ali = countAlignedBasesInRegion(chromAri, chain->tEnd, next->tStart);
		    if(ali >= minBreakingAli)
		    {
			verbose (2, "Non-colinearity at %s:%d-%d because tAli=%d\n",
				 chain->tName, chain->tStart, next->tEnd, ali);
			break;
		    }
		}
		/* Break if number of aligned query bases between chains is >= minBreakingAli */
		if(next->qStart - chain->qEnd >= minBreakingAli)
		{
		    chromAri = hashMustFindVal(qGenomeAri, chain->qName);
		    ali = chain->qStrand == '+' ?
			countAlignedBasesInRegion(chromAri, chain->qEnd, next->qStart) :
			countAlignedBasesInRegion(chromAri, chain->qSize - next->qStart, chain->qSize - chain->qEnd);
		    if(ali >= minBreakingAli)
		    {
			verbose (2, "Non-colinearity at %s:%d-%d because qAli=%d\n",
				 chain->tName, chain->tStart, next->tEnd, ali);
			break;
		    }
		}
	    }
	    else if(queryIsFragmented)
	    { /* No, query chrom differs, but we allow that if the flag queryIsFragmented is set */ 
		/* Break if distance on target genome is >= minBreakingAli */
		if(next->tStart - chain->tEnd >= minBreakingAli)
		{
		    chromAri = hashMustFindVal(tGenomeAri, chain->tName);
		    ali = countAlignedBasesInRegion(chromAri, chain->tEnd, next->tStart);
		    if(ali >= minBreakingAli)
		    {
			verbose (2, "Non-colinearity at %s:%d-%d because tAli=%d (multiple query sequences)\n",
				 chain->tName, chain->tStart, next->tEnd, ali);
			break;
		    }
		}
		/* Calc aligned bases after end of chain1 on its query sequence, and break if value is too high */
		chromAri = hashMustFindVal(qGenomeAri, chain->qName);
		ali = chain->qStrand == '+' ?
		    countAlignedBasesInRegion(chromAri, chain->qEnd, chain->qSize) :
		    countAlignedBasesInRegion(chromAri, 0, chain->qSize - chain->qEnd);
		if(ali >= minBreakingAli)
		{
		    verbose (2, "Non-colinearity at %s:%d-%d because qAli1=%d\n",
			     chain->tName, chain->tStart, next->tEnd, ali);
		    break;
		}
		/* Add aligned bases before start of chain2 on its query sequence, and break if value is too high */
		chromAri = hashMustFindVal(qGenomeAri, next->qName);
		ali += next->qStrand == '+' ?
		    countAlignedBasesInRegion(chromAri, 0, next->qStart) :
		    countAlignedBasesInRegion(chromAri, next->qSize - next->qStart, next->qSize);
		if(ali >= minBreakingAli)
		{
		    verbose (2, "Non-colinearity at %s:%d-%d because qAli1+qAli2=%d\n",
			     chain->tName, chain->tStart, next->tEnd, ali);
		    break; 
		}
	    }
	    else 
	    { /* If the flag isn't set, we break here */
		break;
	    }
	    verbose (1, "Colinearity at %s:%d-%d!\n", chain->tName, chain->tStart, next->tEnd);
	    chain = chain->next;
	}
	/* Break the chain list at the point where colinarity is broken */
	chain->next = NULL;
	chain = next;
    }

    verbose(1, "Made %d synteny blocks from %d chains\n", slCount(sbList), chainCount);

    return sbList;
}


struct bed *syntenyBlockToTargetBed(struct syntenyBlock *sb)
/* Convert syntenyBlock to bed for target genome.  */
{
    struct bed *bed;
    struct chain *chain, *firstChain, *lastChain = NULL;
    struct cBlock *cBlock1, *cBlock2;
    int i, nBedBlocks = 0;
    char qStrand;

    firstChain = sb->chainList;
    qStrand = firstChain->qStrand;
    for(chain = firstChain; chain != NULL; chain = chain->next)
    {
	if(qStrand != chain->qStrand)
	    qStrand = '.';
	nBedBlocks++;
	if(minGapShown != 0)
	    for(cBlock1 = chain->blockList; cBlock1->next != NULL; cBlock1 = cBlock1->next)
		if(cBlock1->next->tStart - cBlock1->tEnd >= minGapShown)
		    nBedBlocks++;
	lastChain = chain;
    }

    AllocVar(bed);
    AllocArray(bed->blockSizes, nBedBlocks);
    AllocArray(bed->chromStarts, nBedBlocks);

    bed->chrom = cloneString(firstChain->tName);
    bed->strand[0] = qStrand;
    bed->chromStart = firstChain->tStart;
    bed->thickStart = firstChain->tStart;
    bed->chromEnd = lastChain->tEnd;
    bed->thickEnd = lastChain->tEnd;
    bed->name = writeLongName ? 
	makeLongQueryLocString(firstChain) :
	makeAbbreviatedLocString(firstChain->qName,
				 firstChain->qStrand == '+' ?
				 firstChain->qStart : 
				 firstChain->qSize - firstChain->qEnd);
    bed->blockCount = nBedBlocks;
    for (i=0, chain = firstChain; chain != NULL; chain = chain->next)
    {
	if(minGapShown == 0)
	{
	    assert(i < nBedBlocks);
	    bed->blockSizes[i] = chain->tEnd - chain->tStart;
	    bed->chromStarts[i] = chain->tStart - firstChain->tStart;
	    i++;
	}
	else
	{
	    cBlock1 = chain->blockList;
	    while(cBlock1)
	    {
		assert(i < nBedBlocks);
		cBlock2 = cBlock1;
		while(cBlock2->next && cBlock2->next->tStart - cBlock2->tEnd < minGapShown)
		    cBlock2 = cBlock2->next;
		bed->blockSizes[i] = cBlock2->tEnd - cBlock1->tStart;
		bed->chromStarts[i] = cBlock1->tStart - firstChain->tStart;
		cBlock1 = cBlock2->next;
		i++;
	    }
	}
    }

    /* COLORS:
       If we want to color features by chromosome, we need to copy code
       from hg/hgTracks/hgTracks.c */

    return bed;
}


struct bed *chainListToQueryBed(struct chain *firstChain, int nChains)
{
    struct bed *bed;
    struct chain *chain, *lastChain = NULL;
    struct cBlock *cBlock1, *cBlock2;
    int i, nBedBlocks = 0;
    char qStrand;

    for(i = 0, chain = firstChain; i < nChains; ++i, chain = chain->next)
    {
	assert(chain != NULL);
	nBedBlocks++;
	if(minGapShown != 0)
	    for(cBlock1 = chain->blockList; cBlock1->next != NULL; cBlock1 = cBlock1->next)
		if(cBlock1->next->qStart - cBlock1->qEnd >= minGapShown)
		    nBedBlocks++;
	lastChain = chain;
    }

    AllocVar(bed);
    AllocArray(bed->blockSizes, nBedBlocks);
    AllocArray(bed->chromStarts, nBedBlocks);
    bed->chrom = cloneString(firstChain->qName);
    bed->strand[0] = qStrand = firstChain->qStrand;
    bed->blockCount = nBedBlocks;
    bed->name = writeLongName ? 
	makeLongTargetLocString(firstChain) :
	makeAbbreviatedLocString(firstChain->tName, firstChain->tStart);

    for (i=0, chain = firstChain; i < nBedBlocks; chain = chain->next)
    {
	if(minGapShown == 0)
	{
	    bed->blockSizes[i] = chain->qEnd - chain->qStart;
	    bed->chromStarts[i] = qStrand == '+' ? 
		chain->qStart - firstChain->qStart :
		lastChain->qEnd - chain->qEnd;
	    i++;
	}
	else
	{
	    cBlock1 = chain->blockList;
	    while(cBlock1)
	    {
		assert(i < nBedBlocks);
		cBlock2 = cBlock1;
		while(cBlock2->next && cBlock2->next->qStart - cBlock2->qEnd < minGapShown)
		    cBlock2 = cBlock2->next;
		bed->blockSizes[i] = cBlock2->qEnd - cBlock1->qStart;
		bed->chromStarts[i] = qStrand == '+' ? 
		    cBlock1->qStart - firstChain->qStart :
		    lastChain->qEnd - cBlock2->qEnd;
		cBlock1 = cBlock2->next;
		i++;
	    }
	}
    }
    
    if(qStrand == '+')
    {
	bed->chromStart = bed->thickStart = firstChain->qStart;
	bed->chromEnd = bed->thickEnd = lastChain->qEnd;
    }
    else
    {
	bed->chromStart = bed->thickStart = firstChain->qSize - lastChain->qEnd;
	bed->chromEnd = bed->thickEnd = firstChain->qSize - firstChain->qStart;
	reverseInts(bed->blockSizes, nBedBlocks);
	reverseInts(bed->chromStarts, nBedBlocks);
    }

    return bed; 
}


struct bed *syntenyBlockToQueryBedList(struct syntenyBlock *sb)
/* Convert syntenyBlock to list of beds for query genome.  */
{
    struct bed *bedList = NULL;
    struct chain *chain1, *chain2;
    int i, nChains;

    chain1 = sb->chainList;
    while(chain1)
    {
	nChains = 1;
	for(chain2 = chain1->next; chain2 && sameString(chain1->qName, chain2->qName); chain2 = chain2->next)
	    nChains++;
	slSafeAddHead(&bedList, chainListToQueryBed(chain1, nChains));
	chain1 = chain2;
    }

    return bedList;
}


void outputSyntenyBlocks(struct syntenyBlock *syntenyBlockList,
			 char *tFileName, char *qFileName)
{
    struct syntenyBlock *sb;
    struct chain *chain;
    struct bed *bed, *bedList;
    FILE *tOut = mustOpen(tFileName, "w");
    FILE *qOut = mustOpen(qFileName, "w");

    if(tTrackOptions)
	fprintf(tOut, "track %s\n",tTrackOptions);
    if(qTrackOptions)
	fprintf(qOut, "track %s\n",qTrackOptions);

    for(sb = syntenyBlockList; sb != NULL; sb = sb->next) 
    {
	bed = syntenyBlockToTargetBed(sb);
	bedTabOutN(bed, 12, tOut);
	bedFree(&bed);
	bed = bedList = syntenyBlockToQueryBedList(sb);
	do
	    bedTabOutN(bed, 12, qOut);
	while(bed = bed->next);
	bedFreeList(&bedList);
    }
    fclose(tOut);
    fclose(qOut);
}


void chainNetSynteny(char *chainFilename, /*char *tSizes, char *qSizes, */
		     char *tNetFilename, char *qNetFilename,
		     char *tOutFilename, char *qOutFilename )
/* chainNetSynteny - Make synteny blocks from alignment nets. */
{
    struct lineFile *chainFile;
    struct chain *chain, *chainList = NULL;
    struct hash *tAri, *qAri;
    struct syntenyBlock *syntenyBlockList;

    /* Build indices of aligned regions from net files */
    verbose(1, "Reading aligned regions from %s...\n", tNetFilename);
    tAri = buildAlignedRegionIndex(tNetFilename);
    verbose(1, "Reading aligned regions from %s...\n", qNetFilename);
    qAri = buildAlignedRegionIndex(qNetFilename);
    
    /* Read through chain file and split each chain into synteny blocks */
    verbose(1, "Reading and splitting chains from %s...\n", chainFilename);
    chainFile = lineFileOpen(chainFilename, TRUE);
    while ((chain = chainRead(chainFile)) != NULL)
    {
	breakChainAtFilledGaps(&chain,
			       hashMustFindVal(tAri, chain->tName),
			       hashMustFindVal(qAri, chain->qName),
			       &chainList);
    }

    /* Find and remove overlaps between synteny blocks */
    verbose(1, "Removing overlap between chains...\n");
    removeChainOverlaps(&chainList);

    /* ... */
    verbose(1, "Creating synteny blocks...\n");
    syntenyBlockList = createSyntenyBlocks(&chainList, tAri, qAri);

    /* Output synteny blocks */
    verbose(1, "Writing result to output files...\n");
    outputSyntenyBlocks(syntenyBlockList, tOutFilename, qOutFilename);

    verbose(1, "Done.\n");
}


int main(int argc, char *argv[])
/* Process command line. */
{
    optionInit(&argc, argv, optionSpecs);
    if (argc != 6)
	usage();
    minAli = optionInt("minAli", minAli);
    minBreakingAli = optionInt("breakingAli", minBreakingAli);
    minGapShown = optionInt("minGapShown", minGapShown);
    queryIsFragmented = optionExists("frag");
    writeLongName = optionExists("longName");
    tTrackOptions = optionVal("tTrackOpt", tTrackOptions);
    qTrackOptions = optionVal("qTrackOpt", qTrackOptions);
    chainNetSynteny(argv[1], argv[2], argv[3], argv[4], argv[5]);

    return 0;
}
