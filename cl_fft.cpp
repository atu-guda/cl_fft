#include <unistd.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <climits>
// #include <complex>
#include <vector>
#include <algorithm>
#include <cerrno>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fftw3.h>

using namespace std;

// using Compl = complex<double>;

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
   explicit Fftw_Data( unsigned a_n_in );
   Fftw_Data( const std::vector<double> &d_x, const std::vector<double> &d_y );
   Fftw_Data( const Fftw_Data &r ) = delete;
   ~Fftw_Data();
   Fftw_Data& operator=( const Fftw_Data &rhs ) = delete;
   double  operator[]( unsigned i, unsigned j ) const { return d[i][j]; };
   double& operator[]( unsigned i, unsigned j )       { return d[i][j]; };
   fftw_complex* data() { return d; };
   const fftw_complex* cdata() { return d; };
   unsigned size() const { return sz; }
   unsigned get_N_in() const { return n_in; }

  protected:
   unsigned n_in, sz;
   fftw_complex *d;
};

Fftw_Data::Fftw_Data( unsigned a_n_in )
  :   n_in( a_n_in ), sz( n_in/2+1 ),
  d( (fftw_complex*) fftw_malloc( sz * sizeof( fftw_complex ) ) )
{
}

Fftw_Data::Fftw_Data( const std::vector<double> &d_x, const std::vector<double> &d_y )
  :   n_in( d_x.size()*2 ), sz( d_x.size() ),
  d( (fftw_complex*) fftw_malloc( (sz+2) * sizeof( fftw_complex ) ) )
{
  for( size_t i=0; i<sz; ++i ) {
    d[i][0] = d_x[i];
    d[i][1] = d_y[i];
  }
}

Fftw_Data::~Fftw_Data()
{
  fftw_free( d );
}

struct pgm_params {
  double f_max     = DBL_MAX;
  double dt        = 0; // 0 means auto
  bool out_complex = false;
  bool in_complex  = false;
  bool drop_zero   = false;
  bool out_Hz      = false;
  int  f_dir       = FFTW_FORWARD;
  int t_idx        =   0;
  int re_idx       =   1;
  int im_idx       =  -1;
  unsigned o_w     =  18; // output width
};

void out_res( ostream &os, const Fftw_Data &d, const pgm_params &p );
size_t read_infile( istream &is, pgm_params &prm, vector<double> &d_x, vector<double> &d_y );
void show_help( const char *pname );

int main( int argc, char **argv )
{
  pgm_params p_prm;
  const char *ofile = nullptr;

  ostream *os = &cout;

  int op;
  while( (op = getopt( argc, argv, "hdcCr0zt:x:f:o:w:T:") ) != -1 ) {
    switch ( op ) {
      case 'h': show_help( argv[0] );           return 0;
      case 'd': debug++;                        break;
      case 'c': p_prm.out_complex = true;       break;
      case 'C': p_prm.in_complex  = true;       break;
      case 'r': p_prm.f_dir = FFTW_BACKWARD;    break;
      case 't': p_prm.t_idx  = atoi( optarg );  break;
      case 'x': p_prm.re_idx = atoi( optarg );  break;
      case 'y': p_prm.im_idx = atoi( optarg );  break;
      case 'w': p_prm.o_w    = atoi( optarg );  break;
      case 'f': p_prm.f_max  = atof( optarg );  break;
      case 'T': p_prm.dt     = atof( optarg );  break;
      case 'o': ofile = optarg;                 break;
      case '0': p_prm.drop_zero = true;         break;
      case 'z': p_prm.out_Hz = true;            break;
      default: cerr << "Unknown or bad option '" << (char)optopt << endl;
               show_help( argv[0] );
               return 1;
    };
  };

  if( p_prm.in_complex && p_prm.im_idx < 0 ) {
    p_prm.im_idx = p_prm.re_idx + 1;
  }

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

  vector<double> in_x, in_y;

  auto n =  read_infile( ifs, p_prm, in_x, in_y );
  if( n < 1 ) {
    cerr << "input data error" << endl;
    return 4;
  }

  ifs.close();

  if( p_prm.f_dir == FFTW_BACKWARD ) { // TODO: implement
    cerr << "backword transform is unimplemented for now" << endl;
    return 10;
    // Fftw_Data out( n, p_prm.d_t );
    // plan =  fftw_plan_r2r_1d( n, &(in_x[0]), out.data(), FFTW_REDFT00, FFTW_ESTIMATE );
    // fftw_execute( plan );

    // out_res( *os, out, p_prm );
  }

  Fftw_Data out( p_prm.in_complex ? 2*n : n );

  if( p_prm.in_complex ) {
    if( in_x.size() != in_y.size() ) {
      cerr << "Error: Re and Im array sizes not equal: " << in_x.size() << " != " << in_y.size() << endl;
      return 7;
    }
    Fftw_Data in_d( in_x, in_y );
    cerr << "# TODO:" << endl;
    fftw_plan plan = fftw_plan_dft_1d( n, in_d.data(), out.data(), FFTW_FORWARD, FFTW_ESTIMATE );
    cerr << "# plan:" << endl;
    fftw_execute( plan );
    cerr << "# exec:" << endl;
    fftw_destroy_plan( plan );
    out_res( *os, out, p_prm );
    return 0;
  }

  // real input
  fftw_plan plan = fftw_plan_dft_r2c_1d( n, &(in_x[0]), out.data(), FFTW_ESTIMATE );

  fftw_execute( plan );
  fftw_destroy_plan( plan );

  out_res( *os, out, p_prm );

  return 0;
}

size_t read_infile( istream &is, pgm_params &prm, vector<double> &d_x, vector<double> &d_y )
{
  size_t n {0}, n_line {0};
  double old_t { 0 }, t_0 {0};

  const auto idx_max = max( max( prm.t_idx, prm.re_idx) , prm.im_idx );
  const size_t init_sz = 1024 * 128; // large enough
  const bool calc_dt = prm.dt <= 0;

  d_x.clear(); d_y.clear();
  d_x.reserve( init_sz );
  if( prm.in_complex ) {
    d_y.reserve( init_sz );
  }

  while( is ) {
    string s;
    getline( is, s );
    s = trim( s );
    ++n_line;
    if( s.empty() || s[0] == '#' || s[0] == ';' ) {
      continue;
    }

    istringstream is( s );

    double c_t {0}, x {0}, y {0};
    for( decltype(+idx_max) i=0; i<= idx_max ; ++i ) {
      string sv;
      is >> sv;
      if( !is || sv.empty() ) {
        cerr << "Read only " << i << " columns in line " << n_line << " n= " << n << endl;
        cerr << "Line: \"" << s << "\"" << endl;
        return 0;
      }

      if( i != prm.t_idx  &&  i != prm.re_idx &&  i != prm.im_idx  ) { // skip unused, TODO: complex
        continue;
      }

      size_t epos;
      double v;
      try {
        v = stod( sv, &epos );
      } catch( std::invalid_argument &e ) {
        cerr << "Fail to convert column " << i << " in line " << n_line << " n= " << n << endl;
        cerr << "Column: \"" << sv << "\" Line: \"" << s << "\"" << e.what() << endl;
        return 0;
      }
      if( epos < 1 ) { // may be overkill
        cerr << "Fail to convert column " << i << " in line " << n_line << " n= " << n << endl;
        return 0;
      }

      if( i == prm.t_idx ) {
        c_t = v;
      } else if( i == prm.re_idx ) {
        x = v;
      } else if( i == prm.im_idx ) {
        y = v;
      } else {
        cerr << "# warn: unused index " << i << endl;
      }
    }

    if( n == 0 ) {
      t_0 = c_t;
    }
    auto c_dt = c_t - old_t;

    if( calc_dt ) {
      if( n == 1 ) {
        if( c_dt <= 0 ) {
          cerr << "Error: Bad delta t value: " << c_dt << " = " << c_t << " - " << old_t;
          return 0;
        }
        prm.dt = c_dt;
      } else {
        if( fabs( ( c_dt - prm.dt ) / prm.dt ) > 3e-1 ) {
          cerr << "Inconsistent dt " << c_dt << ' ' << prm.dt << " in line " << n_line << " n= " << n << endl;
        }
      }
    }
    old_t = c_t;

    d_x.push_back( x );
    if( prm.in_complex ) {
      d_y.push_back( y );
    }
    ++n;
  };

  cerr << "# debug: dt_1= " << prm.dt << " dt_n = " << ( ( old_t - t_0 ) / (n-1) )
       << " sz: " << d_x.size() << ' ' << n << endl;

  return d_x.size();
}

void out_res( ostream &os, const Fftw_Data &d, const pgm_params &p )
{
  const auto o_n = d.size();
  unsigned st = p.drop_zero ? 1 : 0;
  const auto n = d.get_N_in();

  double f_coeff = ( p.out_Hz ? 1 : (2 * M_PI) )  /  ( p.dt * n );

  for( decltype(+o_n) i=st; i<o_n ; ++i ) {
    double fr = f_coeff * i;
    if( fr > p.f_max ) {
      break;
    }
    os << setw( p.o_w) << fr << ' ';
    if( p.out_complex ) {
       os << setw( p.o_w ) << d[i,0] << ' ' << setw( p.o_w) << d[i,1] << endl;
    } else {
       os << setw( p.o_w ) << ( 2.0 * hypot( d[i,0], d[i,1] ) / ( double(n) ) ) << endl;
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
          "  -w width = output number width,def = 18 \n"
          "  -t idx  = index of 't' column, def = 0\n"
          "  -x idx  = index of 'x' column, def = 1\n"
          "  -y idx  = index of 'img' column, def = x+1\n"
          "  -f val  = maximum required frequency \n"
          "  -T dt   = force time step \n"
          "  -o file  = set output file \n"
          "  -0      = drop zero frequency from output \n"
          "  -z      = output in ordinary frequency (Hz), not in omega\n"
          ;
}

