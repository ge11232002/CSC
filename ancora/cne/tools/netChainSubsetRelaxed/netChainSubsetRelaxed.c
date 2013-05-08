/* netChainSubsetRelaxed - Create chain file with subset of chains that appear in the net.
 * Adapted from netChainSubset in the UCSC source.
 * The only difference compared to the original UCSC code is that this program
 * will not abort if a chain listed in the net file isn't found in the chain
 * file. Instead the program will simply skip that net and continue.
 */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "chain.h"
#include "chainNet.h"

char *type = NULL;
boolean splitOnInsert = FALSE;
boolean wholeChains = FALSE;

void usage()
/* Explain usage and exit. */
{
errAbort(
  "netChainSubsetRelaxed - Create chain file with subset of chains that appear in the net\n"
  "usage:\n"
  "   netChainSubsetRelaxed in.net in.chain out.chain\n"
  "options:\n"
  "   -gapOut=gap.tab - Output gap sizes to file\n"
  "   -type=XXX - Restrict output to particular type in net file\n"
  "   -splitOnInsert - Split chain when get an insertion of another chain\n"
  "   -wholeChains - Write entire chain references by net, don't split\n"
  "    when a high-level net is encoundered.  This is useful when nets\n"
  "    have been filtered.\n"
  );
}

struct optionSpec options[] = {
   {"gapOut", OPTION_STRING},
   {"type", OPTION_STRING},
   {"splitOnInsert", OPTION_BOOLEAN},
   {"wholeChains", OPTION_BOOLEAN},
   {NULL, 0},
};


struct chain *chainLookupIfExists(struct hash *hash, int id)
/* Find chain in hash. Return NULL if not found. */
{
char nameBuf[16];
safef(nameBuf, sizeof(nameBuf), "%x", id);
return hashFindVal(hash, nameBuf);
}


void gapWrite(struct chain *chain, FILE *f)
/* Write gaps to simple two column file. */
{
struct cBlock *a, *b;
a = chain->blockList;
for (b = a->next; b != NULL; b = b->next)
    {
    fprintf(f, "%d\t%d\n", b->tStart - a->tEnd, b->qStart - a->qEnd);
    a = b;
    }
}

void writeChainPart(struct chain *chain, int tStart, int tEnd, FILE *f,
	FILE *gapFile)
/* Write out part of a chain. */
{
struct chain *subChain, *chainToFree;

chainSubsetOnT(chain, tStart, tEnd, &subChain, &chainToFree);
if(subChain != NULL)
   {
   chainWrite(subChain, f);
   if (gapFile != NULL)
       gapWrite(subChain, gapFile);
   chainFree(&chainToFree);
   }
}

void writeChainWhole(struct chain *chain, FILE *f, FILE *gapFile)
/* Write out entire chain. */
{
chainWrite(chain, f);
if (gapFile != NULL)
    gapWrite(chain, gapFile);
}

struct cnFill *nextGapWithInsert(struct cnFill *gapList)
/* Find next in list that has a non-empty child.   */
{
struct cnFill *gap;
for (gap = gapList; gap != NULL; gap = gap->next)
    {
    if (gap->children != NULL)
        break;
    }
return gap;
}

void splitWrite(struct cnFill *fill, struct chain *chain, 
    FILE *f, FILE *gapFile)
/* Split chain into pieces if it has inserts.  Write out
 * each piece. */
{
int tStart = fill->tStart;
struct cnFill *child = fill->children;

for (;;)
    {
    child = nextGapWithInsert(child);
    if (child == NULL)
        break;
    writeChainPart(chain, tStart, child->tStart, f, gapFile);
    tStart = child->tStart + child->tSize;
    child = child->next;
    }
writeChainPart(chain, tStart, fill->tStart + fill->tSize, f, gapFile);
}

void convertFill(struct cnFill *fill, 
	struct chain *chain, FILE *f, FILE *gapFile)
/* Convert subset of chain as defined by fill to axt. */
{
if (type != NULL)
    {
    if (!sameString(type, fill->type))
        return;
    }
if (splitOnInsert)
    splitWrite(fill, chain, f, gapFile);
else if (wholeChains)
    writeChainWhole(chain, f, gapFile);
else
    writeChainPart(chain, fill->tStart, fill->tStart + fill->tSize, f, gapFile);
}


void rConvert(struct cnFill *fillList, 
	struct hash *chainHash, FILE *f, FILE *gapFile)
/* Recursively output chains in net as axt. */
{
struct cnFill *fill;
struct chain *chain;
for (fill = fillList; fill != NULL; fill = fill->next)
    {
    if (fill->chainId)
        {
	chain = chainLookupIfExists(chainHash, fill->chainId);
        if(chain != NULL)
	    convertFill(fill, chain, f, gapFile);
        }
    if (fill->children)
        rConvert(fill->children, chainHash, f, gapFile);
    }
}

void netChainSubset(char *netIn, char *chainIn, char *chainOut)
/* netChainSubset - Create chain file with subset of *
 * chains that appear in the net. */
{
struct hash *chainHash;
struct chainNet *net;
struct lineFile *lf = lineFileOpen(netIn, TRUE);
FILE *f = mustOpen(chainOut, "w");
char *gapFileName = optionVal("gapOut", NULL);
FILE *gapFile = NULL;

if (gapFileName)
    gapFile = mustOpen(gapFileName, "w");
chainHash = chainReadAllWithMeta(chainIn, f);
while ((net = chainNetRead(lf)) != NULL)
    {
    verbose(1, "Processing %s\n", net->name);
    rConvert(net->fillList, chainHash, f, gapFile);
    chainNetFree(&net);
    }
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
type = optionVal("type", type);
splitOnInsert = optionExists("splitOnInsert");
wholeChains = optionExists("wholeChains");
if (argc != 4)
    usage();
netChainSubset(argv[1], argv[2], argv[3]);
return 0;
}
