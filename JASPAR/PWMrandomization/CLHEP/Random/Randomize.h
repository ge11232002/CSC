// $Id: Randomize.h,v 1.3 2003/10/23 21:29:51 garren Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// This file must be included to make use of the HEP Random module
// On some compilers the static instance of the HepRandom generator
// needs to be created explicitly in the client code. The static
// generator is assured to be correctly initialized by including this
// header in the client code.

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
// Gabriele Cosmo - Last change: 13th February 1996
// Ken Smith      - Added Ranshi and DualRand engines: 4th June 1998
//                - Added Ranlux64 and MTwist engines: 14th July 1998
//                - Added Hurd160, Hurd288m and TripleRand 6th Aug 1998
// =======================================================================

#ifndef Rndmze_h
#define Rndmze_h 1

// Including Engines ...

#include "defs.h"
#include "DRand48Engine.h"
#include "DualRand.h"
#include "Hurd160Engine.h"
#include "Hurd288Engine.h"
#include "JamesRandom.h"
#include "MTwistEngine.h"
#include "RandEngine.h"
#include "RanecuEngine.h"
#include "RanluxEngine.h"
#include "Ranlux64Engine.h"
#include "RanshiEngine.h"
#include "TripleRand.h"

// Including distributions ...

//#include "CLHEP/Random/RandBinomial.h"
//#include "CLHEP/Random/RandBreitWigner.h"
//#include "CLHEP/Random/RandChiSquare.h"
//#include "CLHEP/Random/RandExponential.h"
#include "RandFlat.h"
//#include "CLHEP/Random/RandBit.h"
#include "RandGamma.h"
#include "RandGauss.h"
//#include "CLHEP/Random/RandGaussQ.h"
//#include "CLHEP/Random/RandGaussT.h"
//#include "CLHEP/Random/RandGeneral.h"
//#include "CLHEP/Random/RandLandau.h"
//#include "CLHEP/Random/RandPoissonQ.h"
//#include "CLHEP/Random/RandPoissonT.h"
//#include "CLHEP/Random/RandStudentT.h"

namespace CLHEP {

#define HepUniformRand() HepRandom::getTheEngine()->flat()

// On some compilers the static instance of the HepRandom generator
// needs to be created explicitly in the client code (i.e. here).

static int HepRandomGenActive = HepRandom::createInstance();

}  // namespace CLHEP

#ifdef ENABLE_BACKWARDS_COMPATIBILITY
//  backwards compatibility will be enabled ONLY in CLHEP 1.9
using namespace CLHEP;
#endif

#endif
