/**
 * @file profileToAlignmentProfile.cpp
 *
 * @author  Ted Holzman <tholzman@scharp.org>
 *  modified from the Align.cpp file by Paul Edlefsen <pedlefsen@gmail.com>
 * @version 0.1
 * @date 02/01/2012
 * @section LICENSE
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * \brief Create one (or more) "Alignment Profiles" from a fasta file and a
 * a profile
 *
 * @section DESCRIPTION
 * Part of the \a profuse library.
 *
 * Converts a "standard" galosh profile to an alignment profile containing
 * positional information about insertions and deletions.  It gleans this
 * information from a fasta file containing several sequences.  profileToAlignmentProfile
 * aligns these sequences, determines the indels, and combines the information
 * with the given profile to produce an Alignment Profile.
 *
 * The following options are understood:
 *<p><pre>
 * Usage:  profileToAlignmentProfile [options] <profile-file-name> <fasta-file-name>
 *
 * profileToAlignmentProfile options:
 *
 * -h [ --help ]                 output this help message
 * --verbosity arg (=0)          set verbosity level (0 - 3)
 * -i [ --individual ]           output individual alignment profiles instead of
 *                               combined
 * -n [ --nseq ] arg             number of sequences to use (default is ALL)
 * -v [ --viterbi ]              use viterbi algorithm
 * </pre>
 *
 */

#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>    
#include <boost/lexical_cast.hpp>

#include "GenAlignmentProfiles.hpp"

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace galosh;

/**
 * \fn int main(int const argc, char ** argv)
 * \brief main driver.  Parses command line and calls gen_alignment_profiles
 * prints alignment profile output.
 */
int
main ( int const argc, char ** argv )
{
  typedef doublerealspace ProbabilityType;
  typedef bfloat ScoreType;
  typedef bfloat MatrixValueType;

#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  try {
    string config_file;
    
    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic( "Generic options" );
    generic.add_options()
      ( "version,v", "print version string" )
      ( "help,h", "produce help message" )
      ( "config,c", po::value<string>( &config_file )->default_value( "profuse.cfg" ),
        "name of a file of a configuration." )
      ;

    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config( "Configuration" );
    config.add_options()
      ( "profile,p", 
        po::value<string>(),
        "filename: where to find the input profile" )
      ( "alignment_profiles_prefix,o", 
        po::value<string>()->default_value( "profillic_profileToAlignmentProfile" ),
        "filename prefix: where to put the output alignment profiles" )
      ( "fasta,f", 
        po::value<string>(),
        "input sequences, in (unaligned) Fasta format" )
      ("individual,i",
       "output individual alignment profiles instead of average")
      ("individual-filename-suffix-pattern,s",
       po::value<string>()->default_value(  "_$n.aprof" ),
       "pattern for filenames by which to differentiate output individual alignment profiles, in which $n will be replaced by the sequence name and $d with be replaced by the sequence number (default: '_$n.aprof'")
      ("nseq,n",
       po::value<int>(),
       "number of sequences to use (default is ALL)")
      ("viterbi,v", // todo: remove this.  it's just for debugging.
       "use viterbi algorithm")
      ;


    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    // TODO: CHANGE THESE
    //po::options_description lengthadjust_opts( "Lengthadjust options" );
    //lengthadjust_opts.add_options()
    //  ( "proposeInsertingOccupancyThreshold,dms.occupancy_threshold",
    //    po::value<double>(),
    //    "DMS threshold for fraction of sequences inserting/deleting to trigger a model edit of a given position" )
    //  ;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters params;
    
    po::options_description cmdline_options;
    cmdline_options.add( generic ).add( params.m_galosh_options_description ).add( config );//.add( lengthadjust_opts );

    po::options_description config_file_options;
    config_file_options.add( params.m_galosh_options_description ).add( config ); //.add( lengthadjust_opts );

    po::options_description visible( "Basic options" );
    visible.add( generic ).add( config );

    po::positional_options_description p;
    p.add( "profile", 1 );
    p.add( "fasta", 1 );
    p.add( "alignment_profiles_prefix", 1 );
        
    store( po::command_line_parser( argc, argv ).options( cmdline_options ).positional( p ).run(), params.m_galosh_options_map );
    notify( params.m_galosh_options_map );

    // TODO: REMOVE
    //cout << params << endl;
    //cout << endl;

#define USAGE() " " << argv[ 0 ] << " [options] <profile file> <gapless fasta sequences file> [<output alignment profiles file prefix>]"

    // Read in the config file.
    if( config_file.length() > 0 ) {
      ifstream ifs( config_file.c_str() );
      if( !ifs ) {
        if(!params.m_galosh_options_map["config"].defaulted()) {         //TAH 3/13 don't choke if config file was defaulted and is missing
           cout << "Can't open the config file named \"" << config_file << "\"\n";
           return 0;
        } 
      } else {
        store( parse_config_file( ifs, config_file_options ), params.m_galosh_options_map );
        notify( params.m_galosh_options_map );
      }
    }

    if( params.m_galosh_options_map.count( "help" ) > 0 ) {
      cout << "Usage: " << USAGE() << endl;
      cout << visible << "\n";
      return 0;
    }

    if( params.m_galosh_options_map.count( "version" ) ) {
      cout << "profileToAlignmentProfile, version 1.0\n";
      return 0;
    }

    //if( params.m_galosh_options_map.count( "debug" ) ) {
    //  cout << "[DEBUGGING]\n";
    //  return 0;
    //}

    // Required options
    if( ( params.m_galosh_options_map.count( "profile" ) == 0 ) || ( params.m_galosh_options_map.count( "fasta" ) == 0 ) ) {
      cout << "Usage: " << USAGE() << endl;
      return 1;
    }
    
    //if( params.m_galosh_options_map.count( "include-path" ) ) {
    //  cout << "Include paths are: " 
    //       << params.m_galosh_options_map["include-path"].as< vector<string> >() << "\n";
    //}
    //
    //    if (params.m_galosh_options_map.count("input-file"))
    //    {
    //        cout << "Input files are: " 
    //             << params.m_galosh_options_map["input-file"].as< vector<string> >() << "\n";
    //    }

    /// Do the work
    const bool indiv_profiles = ( params.m_galosh_options_map ).count( "individual" ) > 0;
    
    GenAlignmentProfiles<ProbabilityType, ScoreType, MatrixValueType, ResidueType, SequenceResidueType> genAlignProf;
    std::vector<DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> alignment_profiles;
    alignment_profiles = genAlignProf.gen_alignment_profiles( params );

    const std::string output_filename_prefix = ( params.m_galosh_options_map )[ "alignment_profiles_prefix" ].template as<string>();
    const std::string individual_filename_suffix_pattern = ( params.m_galosh_options_map )[ "individual-filename-suffix-pattern" ].template as<string>();

    std::string individual_output_filename = output_filename_prefix + individual_filename_suffix_pattern;
    string const * output_filename_ptr = &individual_output_filename;

    /// Output the results
    if( ( output_filename_ptr == NULL ) ) { //|| ( output_filename.compare( "" ) == 0 ) || ( output_filename.compare( "-" ) == 0 ) ) {
      for( int i = 0; i < alignment_profiles.size(); i++ )
      {
        if( alignment_profiles.size() > 1 ) {
          std:cout << "#" << alignment_profiles[ i ].m_comment << std::endl;
        }
        std::cout << alignment_profiles[ i ];
      }  
    } else {
      for( int i = 0; i < alignment_profiles.size(); i++ )
      {
        if( indiv_profiles ) {
          individual_output_filename = output_filename_prefix + individual_filename_suffix_pattern;
          if( alignment_profiles[ i ].m_comment.length() > 0 ) {
            boost::replace_all( individual_output_filename, "$n", alignment_profiles[ i ].m_comment );
          } else {
            boost::replace_all( individual_output_filename, "$n", boost::lexical_cast<std::string>( i ) );
          }
          boost::replace_all( individual_output_filename, "$d", boost::lexical_cast<std::string>( i ) );
          output_filename_ptr = &individual_output_filename;
        }
      std::ofstream fs ( output_filename_ptr->c_str() );
      if( !fs.is_open() ) {
        // TODO: ?
        cerr << "The alignment profiles output file '" << *output_filename_ptr << "' could not be opened." << endl;
      } else {
          fs << alignment_profiles[ i ];
          fs.close();
          // We print out the output files as a side effect
          cout << *output_filename_ptr << endl;
        }  
      } // End foreach alignment_profile i
    } // End if profile_output_filename_ptr != NULL
    
    return 0; // success
  } catch( std::exception& e ) { /// exceptions thrown by boost stuff
    cerr << "error: " << e.what() << endl;
    return 1;
  } catch( string &err ) {      /// exceptions thrown by GenAlignmentProfiles, etc.
    cerr << "error: " << err << endl;
    return 1;
  } catch( ... ) {               /// anything else
    cerr << "Strange unknown exception" << endl;
    return 1;
  }
    
} // main (..)

