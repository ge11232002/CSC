// $Id: RandEngine.cc,v 1.4.4.4 2005/04/15 16:32:53 garren Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                         --- RandEngine ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - mx is initialised to RAND_MAX: 2nd April 1997
//                - Fixed bug in setSeeds(): 15th September 1997
//                - Private copy constructor and operator=: 26th Feb 1998
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter. Removed RAND_MAX and replaced by
//                  pow(0.5,32.). Flat() returns now 32 bits values
//                  obtained by concatenation: 15th Feb 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// J. Marraffino  - Remove dependence on hepString class  13 May 1999
// M. Fischler    - Rapaired bug that in flat() that relied on rand() to      
//                  deliver 15-bit results.  This bug was causing flat()      
//                  on Linux systems to yield randoms with mean of 5/8(!)     
//                - Modified algorithm such that on systems with 31-bit rand()
//                  it will call rand() only once instead of twice. Feb 2004  
// M. Fischler    - Modified the general-case template for RandEngineBuilder  
//                  such that when RAND_MAX is an unexpected value the routine
//                  will still deliver a sensible flat() random.              
// M. Fischler    - Methods for distrib. instance save/restore  12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//                                                                            
// =======================================================================

#include "defs.h"
#include "RandEngine.h"
#include "Random.h"
#include "engineIDulong.h"
#include <string.h>
#include <cmath>	// for pow()
#include <stdlib.h>	// for int()

namespace CLHEP {

#ifdef NOTDEF 
// The way to test for proper behavior of the RandEngineBuilder
// for arbitrary RAND_MAX, on a system where RAND_MAX is some
// fixed specialized value and rand() behaves accordingly, is 
// to set up a fake RAND_MAX and a fake version of rand() 
// by enabling this block.                               
#undef  RAND_MAX                            
#define RAND_MAX 9382956                    
#include "CLHEP/Random/MTwistEngine.h"      
#include "CLHEP/Random/RandFlat.h"          
MTwistEngine * fakeFlat = new MTwistEngine; 
RandFlat rflat (fakeFlat, 0, RAND_MAX+1);   
int rand() { return (int)rflat(); }         
#endif                                      

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

// number of instances with automatic seed selection
int RandEngine::numEngines = 0;

// Maximum index into the seed table
int RandEngine::maxIndex = 215;

std::string RandEngine::name() const {return "RandEngine";}

RandEngine::RandEngine(long seed) 
: mantissa_bit_32( pow(0.5,32.) )
{
   setSeed(seed,0);
   setSeeds(&theSeed,0);
   seq = 0;
}

RandEngine::RandEngine(int rowIndex, int colIndex)
: mantissa_bit_32( pow(0.5,32.) )
{
   long seeds[2];
   long seed;

   int cycle = abs(int(rowIndex/maxIndex));
   int row = abs(int(rowIndex%maxIndex));
   int col = abs(int(colIndex%2));
   long mask = ((cycle & 0x000007ff) << 20 );
   HepRandom::getTheTableSeeds( seeds, row );
   seed = (seeds[col])^mask;
   setSeed(seed,0);
   setSeeds(&theSeed,0);
   seq = 0;
}

RandEngine::RandEngine()
: mantissa_bit_32( pow(0.5,32.) )
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
   seq = 0;
}

RandEngine::RandEngine(std::istream& is)
: mantissa_bit_32( pow(0.5,32.) )
{
   is >> *this;
}

RandEngine::~RandEngine() {}

RandEngine::RandEngine(const RandEngine &p)
: mantissa_bit_32( pow(0.5,32.) )
{
  // Assignment and copy of RandEngine objects may provoke
  // undesired behavior in a single thread environment.
  
  std::cerr << "!!! WARNING !!! - Illegal operation." << std::endl;
  std::cerr << "- Copy constructor and operator= are NOT allowed on "
	    << "RandEngine objects -" << std::endl;
  *this = p;
}

RandEngine & RandEngine::operator = (const RandEngine &p)
{
  // Assignment and copy of RandEngine objects may provoke
  // undesired behavior in a single thread environment.

  std::cerr << "!!! WARNING !!! - Illegal operation." << std::endl;
  std::cerr << "- Copy constructor and operator= are NOT allowed on "
	    << "RandEngine objects -" << std::endl;
  *this = p;
  return *this;
}

void RandEngine::setSeed(long seed, int)
{
   theSeed = seed;
   srand( int(seed) );
   seq = 0;
}

void RandEngine::setSeeds(const long* seeds, int)
{
  setSeed(seeds ? *seeds : 19780503L, 0);
  theSeeds = seeds;
}

void RandEngine::saveStatus( const char filename[] ) const
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
   if (!outFile.bad()) {
     outFile << theSeed << std::endl;
     outFile << seq << std::endl;
   }
#endif
}

void RandEngine::restoreStatus( const char filename[] )
{
   // The only way to restore the status of RandEngine is to
   // keep track of the number of shooted random sequences, reset
   // the engine and re-shoot them again. The Rand algorithm does
   // not provide any way of getting its internal status.

   std::ifstream inFile( filename, std::ios::in);
   if (!checkFile ( inFile, filename, engineName(), "restoreStatus" )) {
     std::cout << "  -- Engine state remains unchanged\n";	 	  
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
        std::cerr << "\nRandEngine state (vector) description improper."
	       << "\nrestoreStatus has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return;
      }
      v.push_back(xin);
    }
    getState(v);
    return;
  }

   long count;
   
   if (!inFile.bad() && !inFile.eof()) {
//     inFile >> theSeed;  removed -- encompased by possibleKeywordInput
     inFile >> count;
     setSeed(theSeed,0);
     seq = 0;
     while (seq < count) flat();
   }
}

void RandEngine::showStatus() const
{
   std::cout << std::endl;
   std::cout << "---------- Rand engine status ----------" << std::endl;
   std::cout << " Initial seed  = " << theSeed << std::endl;
   std::cout << " Shooted sequences = " << seq << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

// ====================================================
// Implementation of flat() (along with needed helpers)
// ====================================================

// Here we set up such that **at compile time**, the compiler decides based on  
// RAND_MAX how to form a random double with 32 leading random bits, using      
// one or two calls to rand().  Some of the lowest order bits of 32 are allowed 
// to be as weak as mere XORs of some higher bits, but not to be always fixed.  
//                                                                              
// The method decision is made at compile time, rather than using a run-time    
// if on the value of RAND_MAX.  Although it is possible to cope with arbitrary 
// values of RAND_MAX of the form 2**N-1, with the same efficiency, the         
// template techniques needed would look mysterious and perhaps over-stress     
// some older compilers.  We therefore only treat RAND_MAX = 2**15-1 (as on     
// most older systems) and 2**31-1 (as on the later Linux systems) as special   
// and super-efficient cases.  We detect any different value, and use an        
// algorithm which is correct even if RAND_MAX is not one less than a power     
// of 2.  

  template <int> struct RandEngineBuilder {     // RAND_MAX any arbitrary value
  static unsigned int thirtyTwoRandomBits(long& seq) {               
                                                            
  static bool prepared = false;                             
  static unsigned int iT;                                   
  static unsigned int iK;                                   
  static unsigned int iS;                                   
  static int iN;                                            
  static double fS;                                         
  static double fT;                                         
                                                            
  if ( (RAND_MAX >> 31) > 0 )                               
  {                                                         
    // Here, we know that integer arithmetic is 64 bits.    
    if ( !prepared ) {                                      
      iS = RAND_MAX + 1;                     
      iK = 1;                                
//    int StoK = S;                          
      int StoK = iS;                         
      if ( (RAND_MAX >> 32) == 0) {          
        iK = 2;                              
//      StoK = S*S;                          
        StoK = iS*iS;                        
      }                                      
      int m;                                 
      for ( m = 0; m < 64; ++m ) {           
        StoK >>= 1;                          
        if (StoK == 0) break;                
      }                                      
      iT = 1 << m;                           
      prepared = true;                       
    }                                        
    int v = 0;                               
    do {                                     
      for ( int i = 0; i < iK; ++i ) {       
        v = iS*v+rand();  ++seq;                   
      }                                      
    } while (v < iT);                        
    return v & 0xFFFFFFFF;                   
                                             
  }                                          
                                             
  else if ( (RAND_MAX >> 26) == 0 )                                       
  {                                                                       
    // Here, we know it is safe to work in terms of doubles without loss  
    // of precision, but we have no idea how many randoms we will need to 
    // generate 32 bits.                                                  
    if ( !prepared ) {                                                    
      fS = RAND_MAX + 1;                                                  
      double twoTo32 = ldexp(1.0,32);                                     
      double StoK = fS;                                                   
      for ( iK = 1; StoK < twoTo32; StoK *= fS, iK++ ) { }                
      int m;                                                              
      fT = 1.0;                                                           
      for ( m = 0; m < 64; ++m ) {                                        
        StoK *= .5;                                                       
        if (StoK < 1.0) break;                                            
        fT *= 2.0;                                                        
      }                                                                   
      prepared = true;                                                    
    }                                                                     
    double v = 0;                                                         
    do {                                                                  
      for ( int i = 0; i < iK; ++i ) {                                    
        v = fS*v+rand(); ++seq;                                                 
      }                                                                   
    } while (v < fT);                                                     
    return ((unsigned int)v) & 0xFFFFFFFF;                                
                                                                          
  }                                                                       
  else                                                                    
  {                                                                       
    // Here, we know that 16 random bits are available from each of       
    // two random numbers.                                                
    if ( !prepared ) {                                                    
      iS = RAND_MAX + 1;                                                  
      int SshiftN = iS;                                                   
      for (iN = 0; SshiftN > 1; SshiftN >>= 1, iN++) { }                  
      iN -= 17;                                                           
    prepared = true;                                                      
    }                                                                     
    unsigned int x1, x2;                                                  
    do { x1 = rand(); ++seq;} while (x1 < (1<<16) );                            
    do { x2 = rand(); ++seq;} while (x2 < (1<<16) );                            
    return x1 | (x2 << 16);                                               
  }                                                                       
                                                                          
  }                                                                       
};                                                                        
                                                                          
template <> struct RandEngineBuilder<2147483647> { // RAND_MAX = 2**31 - 1
  inline static unsigned int thirtyTwoRandomBits(long& seq) {                      
    unsigned int x = rand() << 1; ++seq; // bits 31-1                      
    x ^= ( (x>>23) ^ (x>>7) ) ^1;        // bit 0 (weakly pseudo-random)   
    return x & 0xFFFFFFFF;               // mask in case int is 64 bits    
    }                                                                     
};                                                                        
                                                                          
                                                                           
template <> struct RandEngineBuilder<32767> { // RAND_MAX = 2**15 - 1      
  inline static unsigned int thirtyTwoRandomBits(long& seq) {                       
    unsigned int x = rand() << 17; ++seq; // bits 31-17                      
    x ^= rand() << 2;              ++seq; // bits 16-2                       
    x ^= ( (x>>23) ^ (x>>7) ) ^3;         // bits  1-0 (weakly pseudo-random)
    return x & 0xFFFFFFFF;                // mask in case int is 64 bits     
    }                                                                      
};                                                                         
                                                                           
double RandEngine::flat()                                      
{                                                              
  double r;                                                    
  do { r = RandEngineBuilder<RAND_MAX>::thirtyTwoRandomBits(seq);
     } while ( r == 0 ); 
  return r/4294967296.0; 
}  

void RandEngine::flatArray(const int size, double* vect)
{
   int i;

   for (i=0; i<size; ++i)
     vect[i]=flat();
}

RandEngine::operator unsigned int() {
  return RandEngineBuilder<RAND_MAX>::thirtyTwoRandomBits(seq);
}

std::ostream & RandEngine::put ( std::ostream& os ) const
{
     char beginMarker[] = "RandEngine-begin";
     char endMarker[]   = "RandEngine-end";

     os << " " << beginMarker << "\n";
     os << theSeed << " " << seq << " ";
     os << endMarker << "\n";
     return os;
}

std::vector<unsigned long> RandEngine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<RandEngine>());
  v.push_back(static_cast<unsigned long>(theSeed));
  v.push_back(static_cast<unsigned long>(seq));
  return v;
}

std::istream & RandEngine::get ( std::istream& is )
{
   // The only way to restore the status of RandEngine is to
   // keep track of the number of shooted random sequences, reset
   // the engine and re-shoot them again. The Rand algorithm does
   // not provide any way of getting its internal status.
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"RandEngine-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cout << "\nInput stream mispositioned or"
	       << "\nRandEngine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
  }
  return getState(is);
}

std::string RandEngine::beginTag ( )  { 
  return "RandEngine-begin"; 
}
  
std::istream & RandEngine::getState ( std::istream& is )
{
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nRandEngine state (vector) description improper."
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
  long count;
  is >> count;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"RandEngine-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nRandEngine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
   }
   setSeed(theSeed,0);
   while (seq < count) flat();
   return is;
}

bool RandEngine::get (const std::vector<unsigned long> & v) {
  if (v[0] != engineIDulong<RandEngine>()) {
    std::cerr << 
    	"\nRandEngine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool RandEngine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nRandEngine get:state vector has wrong length - state unchanged\n";
    return false;
  }
  theSeed   = v[1];
  int count = v[2];
  setSeed(theSeed,0);
  while (seq < count) flat();  
  return true;
}
}  // namespace CLHEP
