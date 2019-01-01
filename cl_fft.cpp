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

class Fftw_Data {
  public:
   Fftw_Data( unsigned asz, unsigned a_n_in, double a_dt );
   Fftw_Data( const Fftw_Data &r ) = delete;
   ~Fftw_Data();
   double  operator()( unsigned i, unsigned j ) const { return d[i][j]; };
   double& operator()( unsigned i, unsigned j )       { return d[i][j]; };
   fftw_complex* data() { return d; };
   const fftw_complex* cdata() { return d; };
   unsigned size() const { return sz; }
   unsigned get_N_in() const { return n_in; }
   double get_dt() const { return dt; }

  protected:
   unsigned sz, n_in;
   double dt;
   fftw_complex *d;
};

Fftw_Data::Fftw_Data( unsigned asz, unsigned a_n_in, double a_dt )
  : sz( asz ),  n_in( a_n_in ), dt( a_dt ),
  d( (fftw_complex*) fftw_malloc( sz * sizeof( fftw_complex ) ) )
{
}

Fftw_Data::~Fftw_Data()
{
  fftw_free( d );
}

struct out_params {
  double f_max = DBL_MAX;
  bool out_complex = false;
  bool drop_zero = false;
  bool out_Hz = false;
};

void out_res( ostream &os, const Fftw_Data &d, const out_params &p );
void show_help( const char *pname );

int main( int argc, char **argv )
{
  int op, n = 0, t_idx = 0, x_idx = 1, idx_max;
  int o_n;
  unsigned f_dir = FFTW_FORWARD;
  out_params o_prm;
  const char *ofile = nullptr;

  ostream *os = &cout;

  double dt = 0, old_t = 0;
  vector<double> in_x;

  while( (op = getopt( argc, argv, "hdcrt:x:f:o:0z") ) != -1 ) {
    switch ( op ) {
      case 'h': show_help( argv[0] );  return 0;
      case 'd': debug++;                break;
      case 'c': o_prm.out_complex = true;       break;
      case 'r': f_dir = FFTW_BACKWARD;  break;
      case 't': t_idx = atoi( optarg ); break;
      case 'x': x_idx = atoi( optarg ); break;
      case 'f': o_prm.f_max = atof( optarg ); break;
      case 'o': ofile = optarg;         break;
      case '0': o_prm.drop_zero = true;       break;
      case 'z': o_prm.out_Hz = true;          break;
      default: cerr << "Unknown or bad option '" << (char)optopt << endl;
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


  const char *ifile = argv[optind];
  if( ifile[0] == '-' || ifile[1] == '\0' ) {
    ifile = "/dev/stdin";
  }
  ifstream ifs( ifile );
  if( ! ifs ) {
    cerr << "Fail to open file <" << ifile << "> : " << strerror(errno) << endl;
    return 2;
  }

  ofstream ofs;
  if( ofile ) {
    ofs.open( ofile );
    if( ! ofs ) {
      cerr << "Fail to open output file <" << ofile << "> : " << strerror(errno) << endl;
      return 3;
    }
    os = &ofs;
  }

  idx_max = max( t_idx, x_idx );
  vector<double> vals( idx_max+2 );

  in_x.reserve( 1024 * 128 ); // large enough
  while( ifs ) {
    string s;
    getline( ifs, s );
    if( s.empty() ) {
      continue;
    }
    if( s[0] == '#' || s[0] == ';' ) {
      continue;
    }
    istringstream is( s );
    vals.assign( vals.size(), 0 );

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
  cerr << "# n= " << n << " o_n= " << o_n << " dt = " << dt << endl;

  Fftw_Data out( o_n, n, dt );

  fftw_plan plan;
  plan = fftw_plan_dft_r2c_1d( n, &(in_x[0]), out.data(), FFTW_ESTIMATE );

  fftw_execute( plan );
  fftw_destroy_plan( plan );

  out_res( *os, out, o_prm );

  return 0;
}

void out_res( ostream &os, const Fftw_Data &d, const out_params &p )
{
  unsigned o_n = d.size();
  unsigned st = p.drop_zero ? 1 : 0;
  auto n = d.get_N_in();
  double f_coeff = ( p.out_Hz ? 1 : (2 * M_PI) ) / ( d.get_dt() * n );
  for( unsigned i=st; i<o_n ; ++i ) {
    double fr = f_coeff * i;
    if( fr > p.f_max ) {
      break;
    }
    if( p.out_complex ) {
       os << fr << ' ' << d( i, 0 ) << ' ' << d( i, 1 ) << endl;
    } else {
       os << fr << ' ' << 2.0 * hypot( d( i, 0 ), d( i, 1 ) ) / ( double(n) ) << endl;
    }
  }
}

void show_help( const char *pname )
{
  cerr << "Usage: " << pname << " [opts] infile|-\n opts:\n";
  cerr << "  -h = help\n";
  cerr << "  -d = debug++\n";
  cerr << "  -c = complex output\n";
  cerr << "  -r = reverse transform (unimplemented)\n";
  cerr << "  -t idx  = index of 't' column, def = 0\n";
  cerr << "  -x idx  = index of 'x' column, def = 1\n";
  cerr << "  -f val  = maximum required frequency \n";
  cerr << "  -o file  = set output file \n";
  cerr << "  -0      = drop zero frequency from output \n";
  cerr << "  -z      = output in ordinary frequency (Hz), not in omega\n";
}

