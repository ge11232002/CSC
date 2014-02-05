// $Id: Hurd160Engine.cc,v 1.4.4.2 2005/04/15 16:32:53 garren Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                           HEP Random
// Ken Smith      - Initial draft started:                    23rd Jul 1998
//                - Added conversion operators:                6th Aug 1998
// J. Marraffino  - Added some explicit casts to deal with
//                  machines where sizeof(int) != sizeof(long)  22 Aug 1998
// M. Fischler    - Modified use of the various exponents of 2
//                  to avoid per-instance space overhead and
//                  correct the rounding procedure              15 Sep 1998
// 		  - Modified various methods to avoid any
//		    possibility of predicting the next number
//		    based on the last several:  We skip one 
//		    32-bit word per cycle.			15 Sep 1998
//		  - modify word[0] in two of the constructors
//		    so that no sequence can ever accidentally
//		    be produced by differnt seeds.		15 Sep 1998
// J. Marraffino  - Remove dependence on hepStrings class       13 May 1999
// M. Fischler    - Put endl at end of a save                   10 Apr 2001
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods put, get for instance save/restore   12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//		    
// =======================================================================

#include "defs.h"
#include "Random.h"
#include "Hurd160Engine.h"
#include "engineIDulong.h"
#include <string.h>
#include <cmath>	// for ldexp()
#include <stdlib.h>	// for abs(int)

using namespace std;

namespace CLHEP {

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

double Hurd160Engine::twoToMinus_32;
double Hurd160Engine::twoToMinus_53;
double Hurd160Engine::nearlyTwoToMinus_54;

std::string Hurd160Engine::name() const {return "Hurd160Engine";}

void Hurd160Engine::powersOfTwo() {
  twoToMinus_32 = ldexp (1.0, -32);
  twoToMinus_53 = ldexp (1.0, -53);
  nearlyTwoToMinus_54 = ldexp (1.0, -54) - ldexp (1.0, -100);
}

// Number of instances with automatic seed selection
int Hurd160Engine::numEngines = 0;

// Maximum index into the seed table
int Hurd160Engine::maxIndex = 215;

static inline unsigned int f160(unsigned int a, unsigned int b, unsigned int c)
{
  return ( ((b<<2) & 0x7c) | ((a<<2) & ~0x7c) | (a>>30) ) ^ ( (c<<1)|(c>>31) );
}

Hurd160Engine::Hurd160Engine() {
  powersOfTwo();
  int cycle    = abs(int(numEngines/maxIndex));
  int curIndex = abs(int(numEngines%maxIndex));
  long mask = ((cycle & 0x007fffff) << 8);
  long seedlist[2];
  HepRandom::getTheTableSeeds( seedlist, curIndex );
  seedlist[0] ^= mask;
  seedlist[1] = 0;
  setSeeds(seedlist, 0);
  words[0] ^= 0x1324abcd;	 // To make unique vs long or two unsigned
  if (words[0]==0) words[0] = 1; // ints in the constructor
  ++numEngines;
  for( int i=0; i < 100; ++i ) flat();            // warm-up just a bit
}

Hurd160Engine::Hurd160Engine( std::istream& is ) {
    is >> *this;
}

Hurd160Engine::Hurd160Engine( long seed ) {
  powersOfTwo();
  long seedlist[2]={seed,0};
  setSeeds(seedlist, 0);
  words[0] ^= 0xa5482134;	 // To make unique vs default two unsigned
  if (words[0]==0) words[0] = 1; // ints in the constructor
  for( int i=0; i < 100; ++i ) flat();            // warm-up just a bit
}

Hurd160Engine::Hurd160Engine( int rowIndex, int colIndex ) {
  powersOfTwo();
  int cycle = abs(int(rowIndex/maxIndex));
  int   row = abs(int(rowIndex%maxIndex));
  int   col = colIndex & 0x1;
  long mask = (( cycle & 0x000007ff ) << 20 );
  long seedlist[2];
  HepRandom::getTheTableSeeds( seedlist, row );
  // NOTE: is that really the seed wanted (PGC) ??
  seedlist[0] = (seedlist[col])^mask;
  seedlist[1]= 0;
  setSeeds(seedlist, 0);
  for( int i=0; i < 100; ++i ) flat();            // warm-up just a bit
}

Hurd160Engine::~Hurd160Engine() { }

Hurd160Engine::Hurd160Engine( const Hurd160Engine & p ) 
{
  *this = p;
}

Hurd160Engine & Hurd160Engine::operator=( const Hurd160Engine & p ) {
  if( this != &p ) {
    wordIndex = p.wordIndex;
    for( int i=0; i < 5; ++i ) {
      words[i] = p.words[i];
    }
  }
  return *this;
}

void Hurd160Engine::advance() {
     register unsigned int W0, W1, W2, W3, W4;

     W0 = words[0];
     W1 = words[1];
     W2 = words[2];
     W3 = words[3];
     W4 = words[4];
     W1 ^= W0; W0 = f160(W4, W3, W0);
     W2 ^= W1; W1 = f160(W0, W4, W1);
     W3 ^= W2; W2 = f160(W1, W0, W2);
     W4 ^= W3; W3 = f160(W2, W1, W3);
     W0 ^= W4; W4 = f160(W3, W2, W4);
     words[0] = W0 & 0xffffffff;
     words[1] = W1 & 0xffffffff;
     words[2] = W2 & 0xffffffff;
     words[3] = W3 & 0xffffffff;
     words[4] = W4 & 0xffffffff;
     wordIndex = 5;

} // advance();

double Hurd160Engine::flat() {

  if( wordIndex <= 2 ) {        // MF 9/15/98:
                                // skip word 0 and use two words per flat
    advance();
  }

  return   words[--wordIndex] * twoToMinus_32 + // most significant part
     (words[--wordIndex]>>11) * twoToMinus_53 + // fill in rest of bits
                    nearlyTwoToMinus_54;        // make sure non-zero
}

void Hurd160Engine::flatArray( const int size, double* vect ) {
    for (int i = 0; i < size; ++i) {
        vect[i] = flat();
    }
}

void Hurd160Engine::setSeed( long seed, int ) {
  words[0] = (unsigned int)seed;
  for (wordIndex = 1; wordIndex < 5; ++wordIndex) {
    words[wordIndex] = 69607 * words[wordIndex-1] + 54329;
  }
}

void Hurd160Engine::setSeeds( const long* seeds, int ) {
  setSeed( *seeds ? *seeds : 32767, 0 );
  theSeeds = seeds;
}
     
void Hurd160Engine::saveStatus( const char filename[] ) const {
  std::ofstream outFile(filename, std::ios::out);
  if( !outFile.bad() ) {
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
    outFile << std::setprecision(20) << theSeed << " ";
    outFile << wordIndex << " ";
    for( int i = 0; i < 5; ++i ) {
      outFile << words[i] << " ";
    }
    outFile << std::endl;
#endif
}

void Hurd160Engine::restoreStatus( const char filename[] ) {
  std::ifstream inFile(filename, std::ios::in);
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
        std::cerr << "\nHurd160Engine state (vector) description improper."
	       << "\nrestoreStatus has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return;
      }
      v.push_back(xin);
    }
    getState(v);
    return;
  }

  if( !inFile.bad() ) {
//     inFile >> theSeed;  removed -- encompased by possibleKeywordInput
    inFile >> wordIndex;
    for( int i = 0; i < 5; ++i ) {
        inFile >> words[i];
    }
  }
}

void Hurd160Engine::showStatus() const {
  int pr = std::cout.precision(20);
  std::cout << std::endl;
  std::cout << "----------- Hurd engine status ----------" << std::endl;
  std::cout << "Initial seed  = " << theSeed   << std::endl;
  std::cout << "Current index = " << wordIndex << std::endl;
  std::cout << "Current words = " << std::endl;
  for( int i = 0; i < 5 ; ++i ) {
    std::cout << "    " << words[i] << std::endl;
  }
  std::cout << "------------------------------------------" << std::endl;
  std::cout.precision(pr);
}

Hurd160Engine::operator float() {
  if( wordIndex <= 1 ) {        // MF 9/15/98:  skip word 0
    advance();
  }
  return words[--wordIndex ] * twoToMinus_32;
}

Hurd160Engine::operator unsigned int() {
  if( wordIndex <= 1 ) {        // MF 9/15/98:  skip word 0
    advance();
  }
  return words[--wordIndex];
}

std::ostream& Hurd160Engine::put(std::ostream& os) const {
  char beginMarker[] = "Hurd160Engine-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
#ifdef REMOVED 
  char endMarker[]   = "Hurd160Engine-end";
  int pr = os.precision(20);
  os << " " << beginMarker << " ";
  os << theSeed  << " ";
  os << wordIndex << " ";
  for (int i = 0; i < 5; ++i) {
    os << words[i]  << "\n";
  }
  os << endMarker   << "\n ";
  os.precision(pr);
  return os;
#endif
}

std::vector<unsigned long> Hurd160Engine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<Hurd160Engine>());
  v.push_back(static_cast<unsigned long>(wordIndex));
  for (int i = 0; i < 5; ++i) {
    v.push_back(static_cast<unsigned long>(words[i]));
  }
  return v;
}


std::istream& Hurd160Engine::get(std::istream& is) {
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"Hurd160Engine-begin")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nInput mispositioned or"
	      << "\nHurd160Engine state description missing or"
	      << "\nwrong engine type found." << std::endl;
    return is;
  }
  return getState(is);
}

std::string Hurd160Engine::beginTag ( )  { 
  return "Hurd160Engine-begin"; 
}

std::istream& Hurd160Engine::getState(std::istream& is) {
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nHurd160Engine state (vector) description improper."
		<< "\ngetState() has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return is;
      }
      v.push_back(uu);
    }
    getState(v);
    return (is);
  }

//  is >> theSeed;  Removed, encompassed by possibleKeywordInput()

  char endMarker   [MarkerLen];
  is >> wordIndex;
  for (int i = 0; i < 5; ++i) {
    is >> words[i];
  }
  is >> std::ws;
  is.width(MarkerLen);
  is >> endMarker;
  if (strcmp(endMarker,"Hurd160Engine-end")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nHurd160Engine state description incomplete."
	      << "\nInput stream is probably mispositioned now." << std::endl;
    return is;
  }
  return is;
}


bool Hurd160Engine::get (const std::vector<unsigned long> & v) {
  if (v[0] != engineIDulong<Hurd160Engine>()) {
    std::cerr << 
    	"\nHurd160Engine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool Hurd160Engine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nHurd160Engine get:state vector has wrong length - state unchanged\n";
    return false;
  }
  wordIndex = v[1];
  for (int i = 0; i < 5; ++i) {
    words[i] = v[i+2];
  }
  return true;
}


}  // namespace CLHEP
