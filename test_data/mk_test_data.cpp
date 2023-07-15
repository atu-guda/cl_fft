#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <getopt.h>

void show_help( const char *prm_name );

using namespace std;

using Compl = complex<double>;

const char *fo1 = "x_re.txt";
const char *fo2 = "x_co.txt";

int debug {0};

int main( int argc, char **argv )
{
  double f1  { 50.0  };
  double f1i { -0.01 };
  double dt = 1.0e-4;
  size_t n = 100000;
  unsigned o_w = 18;

  int op;
  while( (op = getopt( argc, argv, "hdn:f:t:b:" ) ) != -1 ) {
    switch ( op ) {
      case 'h': show_help( argv[0] );  return 0;
      case 'd': debug++;               break;
      case 'n': n    = atol( optarg ); break;
      case 't': dt   = atof( optarg ); break;
      case 'f': f1   = atof( optarg ); break;
      case 'b': f1i  = atof( optarg ); break;
      default: cerr << "Unknown or bad option '"<< (char)optopt << "'" << endl;
               show_help( argv[0] );
               return 1;
    };
  };

  Compl cf1 { f1, f1i };
  const double wf1 { 2 * M_PI *  f1 };
  const Compl wcf1 { 2 * M_PI * cf1 };

  ofstream os1 ( fo1 );
  if( ! os1 ) {
    cerr << "Fail to open output file <" << fo1 << "> : " << strerror(errno) << endl;
    return 2;
  }
  ofstream os2 ( fo2 );
  if( ! os1 ) {
    cerr << "Fail to open output file <" << fo2 << "> : " << strerror(errno) << endl;
    return 2;
  }

  for( size_t i = 0; i < n; ++i ) {
    const double t = i * dt;
    os1 << setw(o_w) << t << ' ';    os2 << setw(o_w) << t << ' ';
    os1 << setw(o_w) << sin( wf1 * t ) << " 0.0 " << setw(o_w) << cos( wf1 * t ); // 0.0 for fake imag
    const auto s2 = sin( wcf1 * t );
    os2 << setw(o_w) << s2.real() << ' ' << setw(o_w) << s2.imag() << ' ' << abs(s2);

    os1 << endl; os2 << endl;
  }




  return 0;
}


void show_help( const char *pname )
{
  cerr << "Usage: " << pname << " [opts] infile\n opts:\n";
  cerr << "  -h = help\n";
  cerr << "  -d = debug++\n";
  cerr << "  -t dt = time step\n";
  cerr << "  -f F = real frequency\n";
  cerr << "  -b B = imag frequency\n";
}
