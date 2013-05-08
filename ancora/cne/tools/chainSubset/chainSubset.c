/* chainSubset - Get subset of chains by intersecting with a given region. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "chain.h"
#include "chainNet.h"


void usage()
/* Explain usage and exit. */
{
errAbort(
  "chainSubset - Get subset of chains by intersecting with a given region\n"
  "usage:\n"
  "   chainSubset in.chain chrom start end\n"
  "options:\n"
  "   -q - Coords are on query genome\n"
  );
}


boolean filterByQuery;

struct optionSpec options[] = {
   {"q", OPTION_BOOLEAN},
   {NULL, 0},
};


void chainSubset(char *chainIn, char *chrom, int minStart, int maxEnd)
/* netChainSubset - Create chain file with subset of *
 * chains that appear in the net. */
{
    struct chain *subChain, *chainToFree;
    struct lineFile *lf = lineFileOpen(chainIn, TRUE);
    struct chain *chain;
    int minQStart, maxQEnd;

    while ((chain = chainRead(lf)) != NULL)
    {
	if(filterByQuery)
	{
	    if (sameString(chrom, chain->qName)) 
	    {
		if(chain->qStrand == '+') 
		{
		    minQStart = minStart;
		    maxQEnd = maxEnd;
		}
		else {
		    minQStart = chain->qSize - maxEnd;
		    maxQEnd = chain->qSize - minStart;
		}
		if(chain->qStart <= maxQEnd && chain->qEnd >= minQStart)
		{
		    chainSubsetOnQ(chain, minQStart, maxQEnd, &subChain, &chainToFree);
		    if(subChain != NULL) 
		    {
			chainWrite(subChain, stdout);
			chainFree(&chainToFree);
		    }
		}
	    }
	}
	else 
	{
	    if (sameString(chrom, chain->tName) &&
		chain->tStart <= maxEnd && chain->tEnd >= minStart)
	    {
		chainSubsetOnT(chain, minStart, maxEnd, &subChain, &chainToFree);
		if(subChain != NULL) 
		{
		    chainWrite(subChain, stdout);
		    chainFree(&chainToFree);
		}
	    }
	}
    }
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
filterByQuery = optionExists("q");
if (argc != 5)
    usage();
chainSubset(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
return 0;
}
