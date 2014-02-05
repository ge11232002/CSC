// $Id: TripleRand.cc,v 1.4.4.2 2005/04/15 16:32:53 garren Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                           Hep Random
//                       --- TripleRand ---
//                   class implementation file
// -----------------------------------------------------------------------
// A canopy pseudo-random number generator. Using the Tausworthe
// exclusive-or shift register, a simple Integer Coungruence generator, and
// the Hurd 288 total bit shift register, all XOR'd with each other, we
// provide an engine that should be a fairly good "mother" generator.
// =======================================================================
// Ken Smith      - Initial draft started:                    23rd Jul 1998
//                - Added conversion operators:                6th Aug 1998
// J. Marraffino  - Added some explicit casts to deal with
//                  machines where sizeof(int) != sizeof(long)  22 Aug 1998
// M. Fischler    - Modified use of the various exponents of 2
//                  to avoid per-instance space overhead and
//                  correct the rounding procedure              15 Sep 1998
//                - modify constructors so that no sequence can 
//		    ever accidentally be produced by differnt 
//		    seeds, and so that exxcept for the default
//		    constructor, the seed fully determines the 
//		    sequence.              			15 Sep 1998
// M. Fischler	  - Eliminated dependency on hepString		13 May 1999
// M. Fischler	  - Eliminated Taus() and Cong() which were	
//		    causing spurions warnings on SUN CC		27 May 1999
// M. Fischler    - Put endl at end of puts of Tausworth and	10 Apr 2001
//		    integerCong
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods put, get for instance save/restore   12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//
// =======================================================================

#include "TripleRand.h"
#include "defs.h"
#include "engineIDulong.h"
#include <string.h>
#include <cmath>	// for ldexp()

using namespace std;

namespace CLHEP {

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

//********************************************************************
//   TripleRand  
//********************************************************************

// Number of instances with automatic seed selection
int TripleRand::numEngines = 0;

double TripleRand::twoToMinus_32;
double TripleRand::twoToMinus_53;
double TripleRand::nearlyTwoToMinus_54;

std::string TripleRand::name() const {return "TripleRand";}

void TripleRand::powersOfTwo() {
  twoToMinus_32 = ldexp (1.0, -32);
  twoToMinus_53 = ldexp (1.0, -53);
  nearlyTwoToMinus_54 = ldexp (1.0, -54) - ldexp (1.0, -100);
}

TripleRand::TripleRand() 
: tausworthe (1234567 + numEngines + 175321),
  integerCong(69607 * tausworthe + 54329, numEngines),
  hurd(19781127 + integerCong)
{  
  powersOfTwo();
  theSeed = 1234567;
  ++numEngines;
}

TripleRand::TripleRand(long seed)
: tausworthe ((unsigned int)seed + 175321),
  integerCong(69607 * tausworthe + 54329, 1313),
  hurd(19781127 + integerCong)
{
  powersOfTwo();
  theSeed = seed;
}

TripleRand::TripleRand(std::istream & is) 
{
  is >> *this;
}

TripleRand::TripleRand(int rowIndex, int colIndex)
: tausworthe (rowIndex + numEngines * colIndex + 175321),
  integerCong(69607 * tausworthe + 54329, 19),
  hurd(19781127 + integerCong)
{
  powersOfTwo();
  theSeed = rowIndex;
}

TripleRand::~TripleRand() { }

TripleRand::TripleRand( const TripleRand & p )
{
  *this = p;
}

TripleRand & TripleRand::operator=( const TripleRand & p ) {
  if( this != &p ) {
     tausworthe  = p.tausworthe;
     integerCong = p.integerCong;
     hurd        = p.hurd;
  }
  return *this;
}

double TripleRand::flat() {
  unsigned int ic ( integerCong );
  unsigned int t  ( tausworthe  );
  unsigned int h  ( hurd        );
  return ( (t  ^ ic ^ h) * twoToMinus_32 +          // most significant part
               (h >> 11) * twoToMinus_53 +          // fill in remaining bits
                     nearlyTwoToMinus_54	    // make sure non-zero
         );
}

void TripleRand::flatArray(const int size, double* vect) {
  for (int i = 0; i < size; ++i) {
    vect[i] = flat();
  }
}

void TripleRand::setSeed(long seed, int) {
  theSeed = seed;
  tausworthe  = Tausworthe((unsigned int)seed + numEngines + 175321);
  integerCong = IntegerCong(69607 * tausworthe + 54329, numEngines);
  hurd        = Hurd288Engine( 19781127 + integerCong );
}

void TripleRand::setSeeds(const long * seeds, int) {
  setSeed(seeds ? *seeds : 1234567, 0);
  theSeeds = seeds;
}

void TripleRand::saveStatus(const char filename[]) const {
  std::ofstream outFile(filename, std::ios::out);
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
    outFile << std::setprecision(20) << theSeed << " ";
    tausworthe.put ( outFile );
    integerCong.put( outFile);
    outFile << ConstHurd() << std::endl;
#endif
}

void TripleRand::restoreStatus(const char filename[]) {
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
        std::cerr << "\nTripleRand state (vector) description improper."
	       << "\nrestoreStatus has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return;
      }
      v.push_back(xin);
    }
    getState(v);
    return;
  }

  if (!inFile.bad()) {
//     inFile >> theSeed;  removed -- encompased by possibleKeywordInput
    tausworthe.get ( inFile );
    integerCong.get( inFile );
    inFile >> Hurd();
  }
}

void TripleRand::showStatus() const {
  std::cout << std::setprecision(20) << std::endl;
  std::cout <<         "-------- TripleRand engine status ---------" 
								<< std::endl;
  std::cout << "Initial seed           = " << theSeed << std::endl;
  std::cout << "Tausworthe generator   = " << std::endl;
  tausworthe.put( std::cout ); 
  std::cout << "IntegerCong generator  = " << std::endl;
  integerCong.put( std::cout );
  std::cout << "Hurd288Engine generator= " << std::endl << ConstHurd();
  std::cout << std::endl << "-----------------------------------------" 
								<< std::endl;
}

TripleRand::operator float() {
  return (float)
    ( ( integerCong ^ tausworthe ^ hurd ) * twoToMinus_32 
					+ nearlyTwoToMinus_54 );
					// make sure non-zero!
}

TripleRand::operator unsigned int() {
  return integerCong ^ tausworthe ^ hurd;
}

Hurd288Engine & TripleRand::Hurd() 	       { return hurd; }

const Hurd288Engine & TripleRand::ConstHurd() const 
						{ return hurd; }

std::ostream & TripleRand::put (std::ostream & os ) const {
  char beginMarker[] = "TripleRand-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
#ifdef REMOVED 
  char endMarker[]   = "TripleRand-end";
  int pr=os.precision(20);
  os << " " << beginMarker << "\n";
  os << theSeed << "\n";
  tausworthe.put( os );
  integerCong.put( os );
  os << ConstHurd();
  os << " " <<  endMarker  << "\n";
  os.precision(pr);
  return os;
#endif
}

std::vector<unsigned long> TripleRand::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<TripleRand>());
  tausworthe.put(v);
  integerCong.put(v);
  std::vector<unsigned long> vHurd = hurd.put();
  for (unsigned int i = 0; i < vHurd.size(); ++i) {
    v.push_back (vHurd[i]);
  }
  return v;
}

std::istream & TripleRand::get (std::istream & is) {
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"TripleRand-begin")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nInput mispositioned or"
         << "\nTripleRand state description missing or"
         << "\nwrong engine type found." << std::endl;
    return is;
  }
  return getState(is);
}

std::string TripleRand::beginTag ( )  { 
  return "TripleRand-begin"; 
}
  
std::istream & TripleRand::getState (std::istream & is) {
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nTripleRand state (vector) description improper."
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
  tausworthe.get( is );
  integerCong.get( is );
  is >> Hurd();
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"TripleRand-end")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nTripleRand state description incomplete."
         << "\nInput stream is probably mispositioned now." << std::endl;
    return is;
  }
  return is;
}

bool TripleRand::get (const std::vector<unsigned long> & v) {
  if (v[0] != engineIDulong<TripleRand>()) {
    std::cerr << 
    	"\nTripleRand get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  if (v.size() != VECTOR_STATE_SIZE) {
    std::cerr << "\nTripleRand get:state vector has wrong size: " 
    << v.size() << " - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool TripleRand::getState (const std::vector<unsigned long> & v) {
  std::vector<unsigned long>::const_iterator iv = v.begin()+1;
  if (!tausworthe.get(iv)) return false;
  if (!integerCong.get(iv)) return false;
  std::vector<unsigned long> vHurd;
  while (iv != v.end()) {
    vHurd.push_back(*iv++);
  }
  if (!hurd.get(vHurd)) {
    std::cerr << 
    "\nTripleRand get from vector: problem getting the hurd sub-engine state\n";
    return false;
  }     
  return true;
}

//********************************************************************
//   Tausworthe
//********************************************************************

TripleRand::Tausworthe::Tausworthe() {
  words[0] = 1234567;
  for (wordIndex = 1; wordIndex < 4; ++wordIndex) {
    words[wordIndex] = 69607 * words[wordIndex-1] + 54329;
  } 
}

TripleRand::Tausworthe::Tausworthe(unsigned int seed) {
  words[0] = seed;
  for (wordIndex = 1; wordIndex < 4; ++wordIndex) {
    words[wordIndex] = 69607 * words[wordIndex-1] + 54329;
  }
}

TripleRand::Tausworthe::operator unsigned int() {

// Mathematically:  Consider a sequence of bits b[n].  Repeatedly form
// b[0]' = b[127] ^ b[97]; b[n]' = b[n-1].  This sequence will have a very
// long period (2**127-1 according to Tausworthe's work).

// The actual method used relies on the fact that the operations needed to
// form bit 0 for up to 96 iterations never depend on the results of the
// previous ones.  So you can actually compute many bits at once.  In fact
// you can compute 32 at once -- despite 127 - 97 < 32 -- but 24 was used in
// the method used in Canopy, where they wanted only single-precision float
// randoms.  I will do 32 here.

// When you do it this way, this looks disturbingly like the dread lagged XOR
// Fibonacci.  And indeed, it is a lagged Fibonacii, F(4,3, op) with the op
// being the XOR of a combination of shifts of the two numbers.  Although
// Tausworthe asserted excellent properties, I would be scared to death.
// However, the shifting and bit swapping really does randomize this in a
// serious way.

// Statements have been made to the effect that shift register sequences fail
// the parking lot test because they achieve randomness by multiple foldings,
// and this produces a characteristic pattern.  We observe that in this
// specific algorithm, which has a fairly long lever arm, the foldings become
// effectively random.  This is evidenced by the fact that the generator
// does pass the Diehard tests, including the parking lot test.

// To avoid shuffling of variables in memory, you either have to use circular
// pointers (and those give you ifs, which are also costly) or compute at least
// a few iterations at once.  We do the latter.  Although there is a possible
// trade of room for more speed, by computing and saving 256 instead of 128
// bits at once, I will stop at this level of optimization.

// To remind:  Each (32-bit) step takes the XOR of bits [127-96] with bits
// [95-64] and places it in bits [0-31].  But in the first step, we designate
// word0 as bits [0-31], in the second step, word 1 (since the bits it holds
// will no longer be needed), then word 2, then word 3.  After this, the
// stream contains 128 random bits which we will use as 4 valid 32-bit
// random numbers.

// Thus at the start of the first step, word[0] contains the newest (used)
// 32-bit random, and word[3] the oldest.   After four steps, word[0] again
// contains the newest (now unused) random, and word[3] the oldest.
// Bit 0 of word[0] is logically the newest bit, and bit 31 of word[3]
// the oldest.

  if (wordIndex <= 0) {
    for (wordIndex = 0; wordIndex < 4; ++wordIndex) {
      words[wordIndex] = ( (words[(wordIndex+1) & 3] << 1 ) |
                                   (words[wordIndex] >> 31)   )
                       ^ ( (words[(wordIndex+1) & 3] << 31) |
                                   (words[wordIndex] >>  1)   );
    }
  }
  return words[--wordIndex] & 0xffffffff;
}

void TripleRand::Tausworthe::put( std::ostream & os ) const {
  char beginMarker[] = "Tausworthe-begin";
  char endMarker[]   = "Tausworthe-end";

  int pr=os.precision(20);
  os << " " << beginMarker << " ";
  os << std::setprecision(20);
  for (int i = 0; i < 4; ++i) {
    os << words[i] << " ";
  }
  os << wordIndex;
  os << " " <<  endMarker  << " ";
  os << std::endl;
  os.precision(pr);
}

void TripleRand::Tausworthe::put(std::vector<unsigned long> & v) const {
  for (int i = 0; i < 4; ++i) {
    v.push_back(static_cast<unsigned long>(words[i]));
  }
  v.push_back(static_cast<unsigned long>(wordIndex)); 
}

void TripleRand::Tausworthe::get( std::istream & is ) {
  char beginMarker [MarkerLen];
  char endMarker   [MarkerLen];

  is >> std::ws;
  is.width(MarkerLen);  
  is >> beginMarker;
  if (strcmp(beginMarker,"Tausworthe-begin")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nInput mispositioned or"
              << "\nTausworthe state description missing or"
              << "\nwrong engine type found." << std::endl;
  }
  for (int i = 0; i < 4; ++i) {
    is >> words[i];
  }
  is >> wordIndex;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"Tausworthe-end")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nTausworthe state description incomplete."
              << "\nInput stream is probably mispositioned now." << std::endl;
  }
}

bool 
TripleRand::Tausworthe::get(std::vector<unsigned long>::const_iterator & iv){
  for (int i = 0; i < 4; ++i) {
    words[i] = *iv++;
  }
  wordIndex = *iv++;
  return true;
}

//********************************************************************
//   IntegerCong
//********************************************************************

TripleRand::IntegerCong::IntegerCong() 
: state((unsigned int)3758656018U),
  multiplier(66565),
  addend(12341) 
{
}

TripleRand::IntegerCong::IntegerCong(unsigned int seed, int streamNumber)
: state(seed),
  multiplier(65536 + 1024 + 5 + (8 * 1017 * streamNumber)),
  addend(12341)
{
  // As to the multiplier, the following comment was made:
  // We want our multipliers larger than 2^16, and equal to
  // 1 mod 4 (for max. period), but not equal to 1 mod 8 
  // (for max. potency -- the smaller and higher dimension the 
  // stripes, the better)

  // All of these will have fairly long periods.  Depending on the choice
  // of stream number, some of these will be quite decent when considered
  // as independent random engines, while others will be poor.  Thus these
  // should not be used as stand-alone engines; but when combined with a
  // generator known to be good, they allow creation of millions of good
  // independent streams, without fear of two streams accidentally hitting
  // nearby places in the good random sequence.
}

TripleRand::IntegerCong::operator unsigned int() {
  return state = (state * multiplier + addend) & 0xffffffff;
}

void TripleRand::IntegerCong::put( std::ostream & os ) const {
  char beginMarker[] = "IntegerCong-begin";
  char endMarker[]   = "IntegerCong-end";

  int pr=os.precision(20);
  os << " " << beginMarker << " ";
  os << state << " " << multiplier << " " << addend;
  os << " " <<  endMarker  << " ";
  os << std::endl;
  os.precision(pr);
}

void TripleRand::IntegerCong::put(std::vector<unsigned long> & v) const {
  v.push_back(static_cast<unsigned long>(state));
  v.push_back(static_cast<unsigned long>(multiplier));
  v.push_back(static_cast<unsigned long>(addend));
}

void TripleRand::IntegerCong::get( std::istream & is ) {
  char beginMarker [MarkerLen];
  char endMarker   [MarkerLen];

  is >> std::ws;
  is.width(MarkerLen);  
  is >> beginMarker;
  if (strcmp(beginMarker,"IntegerCong-begin")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nInput mispositioned or"
              << "\nIntegerCong state description missing or"
              << "\nwrong engine type found." << std::endl;
  }
  is >> state >> multiplier >> addend;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"IntegerCong-end")) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nIntegerCong state description incomplete."
              << "\nInput stream is probably mispositioned now." << std::endl;
  }
}

bool 
TripleRand::IntegerCong::get(std::vector<unsigned long>::const_iterator & iv) {
  state      = *iv++;
  multiplier = *iv++;
  addend     = *iv++;
  return true;
}

}  // namespace CLHEP
