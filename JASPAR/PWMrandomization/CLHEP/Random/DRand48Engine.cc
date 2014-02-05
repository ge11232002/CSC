// -*- C++ -*-
// $Id: DRand48Engine.cc,v 1.4.4.2 2005/04/15 16:32:53 garren Exp $
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- DRand48Engine ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// G.Cosmo        - Created: 5th September 1995
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - Added srand48(), seed48(), drand48() implementations
//                  for Windows/NT: 6th March 1997
//                - Fixed bug in setSeeds(): 15th September 1997
//                - Private copy constructor and operator=: 26th Feb 1998
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter: 16th Feb 1998
// J.Marraffino   - Remove dependence on hepString class  13 May 1999
// E.Tcherniaev   - More accurate code for drand48() on NT base on
//                  a code extracted from GNU C Library 2.1.3: 8th Nov 2000
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods for distrib. instacne save/restore  12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/8/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//
// =======================================================================

#include "defs.h"
#include "Random.h"
#include "DRand48Engine.h"
#include "RandomFunc.h"
#include "engineIDulong.h"
#include <string.h>

//#define TRACE_IO

namespace CLHEP {

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 
// Number of instances with automatic seed selection
int DRand48Engine::numEngines = 0;

std::string DRand48Engine::name() const {return "DRand48Engine";}

// Maximum index into the seed table
int DRand48Engine::maxIndex = 215;

DRand48Engine::DRand48Engine(long seed)
{
   setSeed(seed,0);
   setSeeds(&theSeed,0);
}

DRand48Engine::DRand48Engine()
{
   long seeds[2];
   long seed;
   int cycle,curIndex;

   cycle = abs(int(numEngines/maxIndex));
   curIndex = abs(int(numEngines%maxIndex));
   numEngines += 1;
   long mask = ((cycle & 0x007fffff) << 8);
   HepRandom::getTheTableSeeds( seeds, curIndex );
   seed = seeds[0]^mask;
   setSeed(seed,0);
   setSeeds(&theSeed,0);
}

DRand48Engine::DRand48Engine(int rowIndex, int colIndex)
{
   long seed;
   long seeds[2];

   int cycle = abs(int(rowIndex/maxIndex));
   int row = abs(int(rowIndex%maxIndex));
   int col = abs(int(colIndex%2));
   long mask = ((cycle & 0x000007ff) << 20);
   HepRandom::getTheTableSeeds( seeds, row );
   seed = (seeds[col])^mask;
   setSeed(seed,0);
   setSeeds(&theSeed,0);
}

DRand48Engine::DRand48Engine(std::istream& is)
{
   is >> *this;
}

DRand48Engine::~DRand48Engine() {}

DRand48Engine::DRand48Engine(const DRand48Engine &p)
{
  // Assignment and copy of DRand48Engine objects may provoke
  // undesired behavior in a single thread environment.
  
  std::cerr << "!!! WARNING !!! - Illegal operation." << std::endl;
  std::cerr << "- Copy constructor and operator= are NOT allowed on "
	    << "DRand48Engine objects -" << std::endl;
  *this = p;
}

DRand48Engine & DRand48Engine::operator = (const DRand48Engine &p)
{
  // Assignment and copy of DRand48Engine objects may provoke
  // undesired behavior in a single thread environment.

  std::cerr << "!!! WARNING !!! - Illegal operation." << std::endl;
  std::cerr << "- Copy constructor and operator= are NOT allowed on "
	    << "DRand48Engine objects -" << std::endl;
  *this = p;
  return *this;
}

void DRand48Engine::setSeed(long seed, int)
{
   srand48( seed );
   theSeed = seed;
}

void DRand48Engine::setSeeds(const long* seeds, int)
{
  setSeed(seeds ? *seeds : 19780503L, 0);
  theSeeds = seeds;
}

void DRand48Engine::saveStatus( const char filename[] ) const
{
   std::ofstream outFile( filename, std::ios::out ) ;

   if (!outFile.bad()) {
    outFile << "Uvec\n";
    std::vector<unsigned long> v = put();
		     #ifdef TRACE_IO
			 std::cout << "Result of v = put() is:\n"; 
		     #endif
    for (unsigned int i=0; i<v.size(); ++i) {
      outFile << v[i] << "\n";
		     #ifdef TRACE_IO
			   std::cout << v[i] << " ";
			   if (i%6==0) std::cout << "\n";
		     #endif
    }
		     #ifdef TRACE_IO
			 std::cout << "\n";
		     #endif
  }

#ifdef REMOVED
   unsigned short dummy[] = { 0, 0, 0 };
   unsigned short* cseed = seed48(dummy);
   if (!outFile.bad()) {
     outFile << theSeed << std::endl;
     for (int i=0; i<3; ++i) {
       outFile << cseed[i] << std::endl;
       dummy[i] = cseed[i];
     }
     seed48(dummy);
   }
#endif
}

void DRand48Engine::restoreStatus( const char filename[] )
{
   std::ifstream inFile( filename, std::ios::in);
   unsigned short cseed[3];

   if (!checkFile ( inFile, filename, engineName(), "restoreStatus" )) {
     std::cerr << "  -- Engine state remains unchanged\n";
     return;
   }
  if ( possibleKeywordInput ( inFile, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long xin;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      inFile >> xin;
	       #ifdef TRACE_IO
	       std::cout << "ivec = " << ivec << "  xin = " << xin << "    ";
	       if (ivec%3 == 0) std::cout << "\n"; 
	       #endif
      if (!inFile) {
        inFile.clear(std::ios::badbit | inFile.rdstate());
        std::cerr << "\nDRand48Engine state (vector) description improper."
	       << "\nrestoreStatus has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return;
      }
      v.push_back(xin);
    }
    getState(v);
    return;
  }

   if (!inFile.bad() && !inFile.eof()) {
     inFile >> theSeed;
     for (int i=0; i<3; ++i)
//     inFile >> theSeed;  removed -- encompased by possibleKeywordInput
     seed48(cseed);
   }
}

void DRand48Engine::showStatus() const
{
   unsigned short dummy[] = { 0, 0, 0 };
   unsigned short* cseed = seed48(dummy);
   std::cout << std::endl;
   std::cout << "-------- DRand48 engine status ---------" << std::endl;
   std::cout << " Initial seed  = " << theSeed << std::endl;
   std::cout << " Current seeds = " << cseed[0] << ", ";
   std::cout                        << cseed[1] << ", ";
   std::cout                        << cseed[2] << std::endl;
   std::cout << "----------------------------------------" << std::endl;
   for (int i=0; i<3; ++i)
     dummy[i] = cseed[i];
   seed48(dummy);
}

double DRand48Engine::flat()
{
   double num = 0.;

   while (num == 0.)
     num = drand48();
   return num;
}

void DRand48Engine::flatArray(const int size, double* vect)
{
   int i;

   for (i=0; i<size; ++i)
     vect[i]=flat();
}

std::ostream & DRand48Engine::put ( std::ostream& os ) const
{
   char beginMarker[] = "DRand48Engine-begin";
   os << beginMarker << "\nUvec\n";
   std::vector<unsigned long> v = put();
   for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
   }
  return os;  

#ifdef REMOVED
   unsigned short dummy[] = { 0, 0, 0 };
   unsigned short* cseed = seed48(dummy);
   char endMarker[] = "DRand48Engine-end";
   os << " " << beginMarker << " ";
   os << theSeed << " ";
   for (int i=0; i<3; ++i) {
     dummy[i] = cseed[i];
     os << cseed[i] << " ";
   }
   os << endMarker << " ";
   seed48(dummy);
   return os;
#endif
}

std::vector<unsigned long> DRand48Engine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<DRand48Engine>());
  unsigned short dummy[] = { 0, 0, 0 };
  unsigned short* cseed = seed48(dummy);
  for (int i=0; i<3; ++i) {
    dummy[i] = cseed[i];
    v.push_back (static_cast<unsigned long>(cseed[i]));
  }
  seed48(dummy);   
  return v;
}

std::istream & DRand48Engine::get ( std::istream& is )
{
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"DRand48Engine-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nDRand48Engine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
  }
  return getState(is);
}

std::string DRand48Engine::beginTag ( )  { 
  return "DRand48Engine-begin"; 
}

std::istream & DRand48Engine::getState ( std::istream& is )  
{
  unsigned short cseed[3];
  if ( possibleKeywordInput ( is, "Uvec", cseed[0] ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    #ifdef TRACE_IO
    std::cout << "DRand48Engine::getState detected Uvec keyword\n";
    uu = 999999;
    #endif
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      uu = 999999;
      is >> uu;
      #ifdef TRACE_IO
      std::cout << "ivec = " << ivec << " uu = " << uu << "\n";
      #endif
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nDRand48Engine state (vector) description improper."
		<< "\ngetState() has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return is;
      }
      v.push_back(uu);
    }
    getState(v);
    return (is);
  }

//  is >> cseed[0] was removed from loop, encompassed by possibleKeywordInput()

  char endMarker   [MarkerLen];
  is >> theSeed;
  for (int i=1; i<3; ++i) {
    is >> cseed[i];
  }
  is >> std::ws;
  is.width(MarkerLen);
  is >> endMarker;
  if (strcmp(endMarker,"DRand48Engine-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nDRand48Engine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
   }
   seed48(cseed);
   return is;
}

bool DRand48Engine::get (const std::vector<unsigned long> & v) {
  if (v[0] != engineIDulong<DRand48Engine>()) {
    std::cerr << 
    	"\nDRand48Engine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool DRand48Engine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
   "\nDRand48Engine getState:state vector has wrong length - state unchanged\n";
    return false;
  }
  unsigned short cseed[3];
  for (int i=0; i<3; ++i) {
    cseed[i] = static_cast<unsigned short>(v[i+1]);
  }
  seed48(cseed);
  return true;
}

}  // namespace CLHEP
