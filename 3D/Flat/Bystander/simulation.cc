#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>

using namespace std;

using boost::mt19937;
using boost::random::uniform_real_distribution;
using boost::random::uniform_int_distribution;

int rownd (double a) {
        return(int(a+ 0.5));
}

// determine winner of cell competition
int winstrain(int a, int b, int c, double sg, double sw, double rannum) {

  double g_a = 1.0 - sw*(double)(a-1)*(double)(a-2)/2.0 - sg*(double)(2-a)*(double)(3-a)/2.0;
  double g_b = 1.0 - sw*(double)(b-1)*(double)(b-2)/2.0 - sg*(double)(2-b)*(double)(3-b)/2.0;
  double g_c = 1.0 - sw*(double)(c-1)*(double)(c-2)/2.0 - sg*(double)(2-c)*(double)(3-c)/2.0;

  vector< vector< double > > vec = {{g_a/(g_a+g_b+g_c),double(a)}, {g_b/(g_a+g_b+g_c),double(b)}, {g_c/(g_a+g_b+g_c),double(c)}};

  double G = 0.0;
  double result;

  for (int s = 0; s < 3; s++) {
    G += vec[s][0];
    if (G > rannum) {
      result = int(vec[s][1]);
      break;
    }
  }

  return result;

}

// calls a%b with the result always positive
int mod(int a, int b)
{
  int r = a%b;
  return r<0? r+b:r;
}

int main(int argc, char* argv[])
{
  int lattsize, nogen;
  double mu, sG, sW;
  double rannum1, rannum2;

  string statsfilename;  
  ofstream outstats;
  ofstream outstats2;
  ofstream outstats3;
  ifstream testoutstats;
  int filecounter=0;
  string tempstr;

  ostringstream tempstring;	

  int seed = time(0);

  mt19937 generator(seed);
  boost::random::uniform_real_distribution< double > dis(0,1);
  boost::random::uniform_int_distribution< int > disi12(1,2);
  boost::random::uniform_int_distribution< int > disi13(1,3);

    lattsize = atoi(argv[1]);
    nogen = atoi(argv[2]);

    tempstring << filecounter;

    statsfilename += "simulation_run";
    statsfilename += tempstring.str();

    testoutstats.open(statsfilename.c_str());
    testoutstats.close();
	
    while (!testoutstats.fail())
    {
      tempstr = tempstring.str();
      statsfilename.erase(statsfilename.end()-tempstr.size(),statsfilename.end());
      filecounter++;
      tempstring.str("");
      tempstring.clear();
      tempstring << filecounter;
      statsfilename += tempstring.str();
      testoutstats.open(statsfilename.c_str());
      testoutstats.close();
    }

    testoutstats.clear(ios::failbit);
    outstats.open(statsfilename.c_str());

    cout << statsfilename.c_str() << endl;

  sG = 0.1;
  sW = 0.3;
  mu = 0.1;

  cout << "# Lx: " << lattsize << " Lt: " << nogen << " sW = " << sW << endl;
  cout << "t\tx\ty\tm" << endl;

  outstats << "# Lx: " << lattsize << " Lt: " << nogen << " sW = " << sW << endl;
  outstats << "t\tx\ty\tm" << endl;

  vector< vector< int > > lattup(lattsize, vector< int >(lattsize));
  vector< vector< int > > lattdown(lattsize, vector< int >(lattsize));
    
    // initialize lattice
    for (int i = 0; i < lattsize; i++) {
      for (int j = 0; j < lattsize; j++) {
        if (mu > 0.0000000001) {
          lattdown[i][j] = disi12(generator);
          outstats << 0 << "\t" << i-0.5*(j%2) << "\t" << j*sqrt(3)/2.0 << "\t" << lattdown[i][j] << endl;
        } else {
          lattdown[i][j] = disi13(generator);
          outstats << 0 << "\t" << i-0.5*(j%2) << "\t" << j*sqrt(3)/2.0 << "\t" << lattdown[i][j] << endl;
        }
      }
    }

    // population growth
    for (int t=1; t<=nogen; t++) {

      // odd step
      for (int i = 0; i < lattsize; i++) {
        for (int j = 0; j < lattsize; j++) {

          if ( j % 2 == 0) {
            if ((lattdown[i][j] == lattdown[i][mod(j+1,lattsize)]) && (lattdown[i][j] == lattdown[mod(i+1,lattsize)][mod(j+1,lattsize)])) {
              lattup[i][j] = lattdown[i][j];
            } else {
              rannum1 = dis(generator);
              lattup[i][j] = winstrain(lattdown[i][j], lattdown[i][mod(j+1,lattsize)], lattdown[mod(i+1,lattsize)][mod(j+1,lattsize)], sG, sW, rannum1);
            }
          } else {
            if ((lattdown[i][j] == lattdown[mod(i-1,lattsize)][mod(j+1,lattsize)]) && (lattdown[i][j] == lattdown[i][mod(j+1,lattsize)])) {
              lattup[i][j] = lattdown[i][j];
            } else {
	      rannum1 = dis(generator);
              lattup[i][j] = winstrain(lattdown[i][j], lattdown[mod(i-1,lattsize)][mod(j+1,lattsize)], lattdown[i][mod(j+1,lattsize)], sG, sW, rannum1);
            }
          }
          if (lattup[i][j] == 2 && mu > 0.0000000001) {
            rannum2 = dis(generator);
            if (rannum2 < mu) lattup[i][j]++;
          }
          outstats << t*sqrt(6)/3.0 << "\t" << i-0.5*(j%2) << "\t" << j*sqrt(3)/2.0+(t%2)/sqrt(3) << "\t" << lattup[i][j] << endl;

        }
      }
      ++t;
  
      // even step
      for (int i = 0; i < lattsize; i++) {
        for (int j = 0; j < lattsize; j++) {

          if ( j % 2 == 0) {
            if ((lattup[i][j] == lattup[i][mod(j-1,lattsize)]) && (lattup[i][j] == lattup[mod(i+1,lattsize)][mod(j-1,lattsize)])) {
              lattdown[i][j] = lattup[i][j];
            } else {
	      rannum1 = dis(generator);
              lattdown[i][j] = winstrain(lattup[i][j], lattup[i][mod(j-1,lattsize)], lattup[mod(i+1,lattsize)][mod(j-1,lattsize)], sG, sW, rannum1);
            }
          } else {
            if ((lattup[i][j] == lattup[i][mod(j-1,lattsize)]) && (lattup[i][j] == lattup[mod(i-1,lattsize)][mod(j-1,lattsize)])) {
              lattdown[i][j] = lattup[i][j];
            } else {
	      rannum1 = dis(generator);
              lattdown[i][j] = winstrain(lattup[i][j], lattup[i][mod(j-1,lattsize)], lattup[mod(i-1,lattsize)][mod(j-1,lattsize)], sG, sW, rannum1);
            }
          }
          if (lattdown[i][j] == 2 && mu > 0.0000000001) {
            rannum2 = dis(generator);
            if (rannum2 < mu) lattdown[i][j]++;
          }
          outstats << t*sqrt(6)/3.0 << "\t" << i-0.5*(j%2) << "\t" << j*sqrt(3)/2.0+(t%2)/sqrt(3) << "\t" << lattdown[i][j] << endl;

        }
      }

    }

  return 0;
  
}
  
  

