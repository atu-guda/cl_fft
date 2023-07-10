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

string trim( const string &s )
{
  auto i = s.begin();
  for( ; i != s.end() && isspace( *i ) ; ++i ) {}
  auto ri = s.rbegin();
  for( ; ri.base() != i && isspace( *ri ); ++ri ) {}
  return string( i, ri.base() );
}

class Fftw_Data {
  public:
   Fftw_Data( unsigned a_n_in, double a_dt );
   Fftw_Data( const Fftw_Data &r ) = delete;
   ~Fftw_Data();
   Fftw_Data& operator=( const Fftw_Data &rhs ) = delete;
   double  operator()( unsigned i, unsigned j ) const { return d[i][j]; };
   double& operator()( unsigned i, unsigned j )       { return d[i][j]; };
   fftw_complex* data() { return d; };
   const fftw_complex* cdata() { return d; };
   unsigned size() const { return sz; }
   unsigned get_N_in() const { return n_in; }
   double get_dt() const { return dt; }

  protected:
   unsigned n_in, sz;
   double dt;
   fftw_complex *d;
};

Fftw_Data::Fftw_Data( unsigned a_n_in, double a_dt )
  :   n_in( a_n_in ), sz( n_in/2+1 ), dt( a_dt ),
  d( (fftw_complex*) fftw_malloc( sz * sizeof( fftw_complex ) ) )
{
}

Fftw_Data::~Fftw_Data()
{
  fftw_free( d );
}

struct pgm_params {
  double f_max     = DBL_MAX;
  bool out_complex = false;
  bool in_complex  = false;
  bool drop_zero   = false;
  bool out_Hz      = false;
};

void out_res( ostream &os, const Fftw_Data &d, const pgm_params &p );
size_t read_infile( istream &is, unsigned t_idx, unsigned x_idx, vector<double> &d, double &dt );
void show_help( const char *pname );

int main( int argc, char **argv )
{
  int op, t_idx = 0, x_idx = 1;
  unsigned f_dir = FFTW_FORWARD;
  pgm_params p_prm;
  const char *ofile = nullptr;

  ostream *os = &cout;

  double dt = 0;
  vector<double> in_x;

  while( (op = getopt( argc, argv, "hdcCrt:x:f:o:0z") ) != -1 ) {
    switch ( op ) {
      case 'h': show_help( argv[0] );           return 0;
      case 'd': debug++;                        break;
      case 'c': p_prm.out_complex = true;       break;
      case 'C': p_prm.in_complex  = true;       break;
      case 'r': f_dir = FFTW_BACKWARD;          break;
      case 't': t_idx = atoi( optarg );         break;
      case 'x': x_idx = atoi( optarg );         break;
      case 'f': p_prm.f_max = atof( optarg );   break;
      case 'o': ofile = optarg;                 break;
      case '0': p_prm.drop_zero = true;         break;
      case 'z': p_prm.out_Hz = true;            break;
      default: cerr << "Unknown or bad option '" << (char)optopt << endl;
               show_help( argv[0] );
               return 1;
    };
  };

  if( optind != argc-1 ) {
    cerr << "Error in parameters: need input filename\n";
    show_help( argv[0] );
    return 1;
  };



  const char *ifile = argv[optind];
  if( ifile[0] == '-' && ifile[1] == '\0' ) {
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

  auto n =  read_infile( ifs, t_idx, x_idx, in_x, dt );
  if( n < 1 ) {
    cerr << "input data error" << endl;
    return 4;
  }

  ifs.close();
  cerr << "# n= " << n << " dt = " << dt << endl;

  if( f_dir == FFTW_BACKWARD ) { // TODO: implement
    cerr << "backword transform is unimplemented for now" << endl;
    return 10;
    // Fftw_Data out( n, dt );
    // plan =  fftw_plan_r2r_1d( n, &(in_x[0]), out.data(), FFTW_REDFT00, FFTW_ESTIMATE );
    // fftw_execute( plan );

    // out_res( *os, out, p_prm );
  }

  if( p_prm.in_complex ) {
    return 11; // TODO
    // return 0;
  }

  // real input
  fftw_plan plan;
  Fftw_Data out( n, dt );

  plan = fftw_plan_dft_r2c_1d( n, &(in_x[0]), out.data(), FFTW_ESTIMATE );

  fftw_execute( plan );
  fftw_destroy_plan( plan );

  out_res( *os, out, p_prm );

  return 0;
}

size_t read_infile( istream &is, unsigned t_idx, unsigned x_idx, vector<double> &d, double &dt )
{
  size_t n {0}, n_line {0};
  double old_t { 0 };

  auto idx_max = max( t_idx, x_idx );

  d.reserve( 1024 * 128 ); // large enough
  while( is ) {
    string s;
    getline( is, s );
    s = trim( s );
    ++n_line;
    if( s.empty() || s[0] == '#' || s[0] == ';' ) {
      continue;
    }

    istringstream is( s );

    double c_t {0}, x {0};
    for( unsigned i=0; i<= idx_max ; ++i ) {
      string sv;
      is >> sv;
      if( !is || sv.empty() ) {
        cerr << "Read only " << i << " columns in line " << n_line << " n= " << n << endl;
        cerr << "Line: \"" << s << "\"" << endl;
        return 0;
      }

      if( i != t_idx  &&  i != x_idx ) { // skip unused, TODO: complex
        continue;
      }

      size_t epos;
      double v = stod( sv, &epos );
      if( epos < 1 ) {
        cerr << "Fail to convert column " << i << " in line " << n_line << " n= " << n << endl;
        cerr << "Column: \"" << sv << "\" Line: \"" << s << "\"" << endl;
        return 0;
      }
      if( i == t_idx ) {
        c_t = v;
      } else if( i == x_idx ) {
        x = v;
      } else {
        cerr << "# warn: unused index " << i << endl;
      }
    }

    auto c_dt = c_t - old_t;
    if( n == 0 ) {
      // NOP
    } else if( n == 1 ) {
      dt = c_dt;
      if( dt <= 0 ) {
        cerr << "Error: Bad delta t value: " << dt << " = " << c_t << " - " << old_t;
        return 0;
      }
    } else {
      if( fabs( ( c_dt - dt ) / dt ) > 3e-1 ) {
        cerr << "Inconsistent dt " << c_dt << ' ' << dt << " in line " << n_line << " n= " << n << endl;
      }
    }
    old_t = c_t;

    d.push_back( x );
    ++n;
  };

  cerr << "# debug: dt_1= " << dt << " dt_n = " << ( old_t / n ) << endl;

  return n;
}

void out_res( ostream &os, const Fftw_Data &d, const pgm_params &p )
{
  const auto o_n = d.size();
  unsigned st = p.drop_zero ? 1 : 0;
  const auto n = d.get_N_in();

  double f_coeff = ( p.out_Hz ? 1 : (2 * M_PI) )  /  ( d.get_dt() * n );

  for( decltype(+o_n) i=st; i<o_n ; ++i ) {
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
  cerr << "Usage: " << pname << " [opts] infile|-\n opts:\n"
          "  -h = help\n"
          "  -d = debug++\n"
          "  -c = complex output\n"
          "  -C = complex input\n"
          "  -r = reverse transform (unimplemented)\n"
          "  -t idx  = index of 't' column, def = 0\n"
          "  -x idx  = index of 'x' column, def = 1, complex: x+1 \n"
          "  -f val  = maximum required frequency \n"
          "  -o file  = set output file \n"
          "  -0      = drop zero frequency from output \n"
          "  -z      = output in ordinary frequency (Hz), not in omega\n"
          ;
}

