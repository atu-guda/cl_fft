#include <unistd.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cerrno>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fftw3.h>

using namespace std;

int debug = 0;
const int max_cols = 16;
void show_help(const char *pname );

int main( int argc, char **argv )
{
  int op, n = 0, t_idx = 0, x_idx = 1, idx_max;
  int o_n;
  unsigned f_dir = FFTW_FORWARD;
  double f_max = DBL_MAX;
  const char *ifile;
  string s;
  bool calc_cmpl = false, drop_zero = false, out_Hz = false;

  double dt = 0, old_t = 0;
  vector<double> in_x;
  fftw_complex *out;
  fftw_plan plan;

  while( (op = getopt( argc, argv, "hdcrt:x:f:0z") ) != -1 ) {
    switch ( op ) {
      case 'h': show_help( argv[0] );  return 0;
      case 'd': debug++;                break;
      case 'c': calc_cmpl = true;       break;
      case 'r': f_dir = FFTW_BACKWARD;  break;
      case 't': t_idx = atoi( optarg ); break;
      case 'x': x_idx = atoi( optarg ); break;
      case 'f': f_max = atof( optarg ); break;
      case '0': drop_zero = true;       break;
      case 'z': out_Hz = true;          break;
      default: fprintf( stderr, "Unknown or bad option <%c>\n", optopt );
               show_help( argv[0] );
               return 1;
    };
  };

  if( optind != argc-1 ) {
    fprintf( stderr, "Error in parameters: need input filename\n" );
    show_help( argv[0] );
    return 1;
  };

  if( f_dir == FFTW_BACKWARD ) { // TODO: implement
    cerr << "backword transform is unimplemented now" << endl;
    return 5;
  }

  ifile = argv[optind];
  if( ifile[0] == '-' || ifile[1] == '\0' ) {
    ifile = "/dev/stdin";
  }
  idx_max = max( t_idx, x_idx );
  vector<double> vals( idx_max+2 );

  ifstream ifs( ifile );
  if( ! ifs ) {
    cerr << "Fail to open file <" << ifile << "> : " << strerror(errno) << endl;
    return 2;
  }

  while( ifs ) {
    getline( ifs, s );
    if( s.empty() ) {
      continue;
    }
    if( s[0] == '#' || s[0] == ';' ) {
      continue;
    }
    istringstream is( s );
    vals.assign( vals.size(), 0 );
    in_x.reserve( 1024 * 128 ); // large enough

    int i; // need after for
    for( i=0; i<= idx_max ;  ) {
      double v = DBL_MAX;
      is >> v;
      if( v == DBL_MAX ) {
        cerr << "Read only " << i << " columns in line " << n << endl;
        return 3;
      }
      vals[i] = v;
      ++i;
    }
    if( n == 0 ) {
      old_t = vals[t_idx];
    }
    if( n == 1 ) {
      dt = vals[t_idx] - old_t;
      if( dt <= 0 ) {
        cerr << "Bad delta t value: " << dt << " = " << vals[t_idx]
             << " - " << old_t;
      }
    }
    in_x.push_back( vals[x_idx] );
    ++n;
  };

  ifs.close();
  o_n = 1 + n/2;
  cerr << "n= " << n << " dt = " << dt << endl;

  out = (fftw_complex*) fftw_malloc( (2+o_n) * sizeof( fftw_complex ) );

  plan = fftw_plan_dft_r2c_1d( n, &(in_x[0]), out, FFTW_ESTIMATE );

  fftw_execute( plan );
  fftw_destroy_plan( plan );

  for( int i=0; i<o_n ; ++i ) {
    if( i == 0 && drop_zero ) {
      continue;
    }
    double fr = (double)(i) / ( dt*n );
    if( !out_Hz ) {
      fr *= 2 * M_PI;
    }
    if( fr > f_max ) {
      break;
    }
    if( calc_cmpl ) {
       cout << fr << ' ' << out[i][0] << ' ' << out[i][1] << endl;
    } else {
       cout << fr << ' ' << 2.0 * hypot( out[i][0], out[i][1] ) / ( double(n) ) << endl;
    }
  }

  return 0;
}

void show_help( const char *pname )
{
  cerr << "Usage: " << pname << " [opts] infile\n opts:\n";
  cerr << "  -h = help\n";
  cerr << "  -d = debug++\n";
  cerr << "  -c = complex output\n";
  cerr << "  -r = reverse transform (unimplemented)\n";
  cerr << "  -t idx  = index of 't' column, def = 0\n";
  cerr << "  -x idx  = index of 'x' column, def = 1\n";
  cerr << "  -f val  = maximum required frequency \n";
  cerr << "  -0      = drop zero frequency from output \n";
  cerr << "  -z      = output in ordinary frequency (Hz), not in omega\n";
}

