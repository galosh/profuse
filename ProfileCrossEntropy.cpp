/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      ProfileCrossEntropy.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The profileCrossEntropy program.  It calculates the crossEntropy
##      between two Profile HMMs (or if you call it with one profile it
##      calculates the self entropy).
##
#******************************************************************************
#*
#*    This file is part of profuse, a suite of programs for working with
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2015, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
#*    profuse is free software: you can redistribute it and/or modify it under
#*    the terms of the GNU Lesser Public License as published by the Free
#*    Software Foundation, either version 3 of the License, or (at your option)
#*    any later version.
#*
#*    profuse is distributed in the hope that it will be useful, but WITHOUT
#*    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#*    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for
#*    more details.
#*
#*    You should have received a copy of the GNU Lesser Public License along
#*    with profuse.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************/

#include "Algebra.hpp"
#include "Profile.hpp"

#include <iostream>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

int
main ( int const argc, char const ** argv )
{
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (galosh Profile) filename 1> [<input (galosh Profile) filename 2>]" << endl;
    exit( 1 );
  }
  const bool be_verbose = ( argc > 3 );
  if( be_verbose ) {
    cout << "Reading profile from file '" << argv[ 1 ] << "'" << endl;
  }
  galosh::ProfileTreeRoot<ResidueType, floatrealspace> profile1;
  profile1.fromFile( argv[ 1 ] );
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << profile1;
    cout << endl;
  }

  if( argc >= 2 ) {
    if( be_verbose ) {
      cout << "Reading another profile from file '" << argv[ 1 ] << "'" << endl;
    }
    galosh::ProfileTreeRoot<ResidueType, floatrealspace> profile2;
    profile2.fromFile( argv[ 2 ] );
    if( be_verbose ) {
      cout << "\tgot:" << std::endl;
      cout << profile2;
      cout << endl;
    }
    cout << "Cross Entropy: " << profile1.crossEntropy( profile2 ) << endl;
  } else {
    cout << "Self Entropy: " << profile1.crossEntropy( profile1 ) << endl;
  }

  exit( 0 );
} // main (..)
