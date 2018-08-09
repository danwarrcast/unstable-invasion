#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>
#include <mpi.h>

using namespace std;

using boost::mt19937;
using boost::uniform_01;

// calls uniform random number between 0 and 1
double random01(mt19937 generator)
{
  static uniform_01<mt19937> dist(generator);
  return dist();
}

int rbetween(mt19937 generator, int min, int max)
{

  if (min == max) return min;

  return floor(random01(generator)*(1+max-min)+min);

}

int rownd (double a) {
        return(int(a+ 0.5));
}

// calls a%b with the result always positive
int mod(int a, int b)
{
  int r = a%b;
  return r<0? r+b:r;
}

double prob(int m1, int m2, double s)
{
  double g1, g2, result;
  if (m1 == 1) {
    g1 = 1;
  } else if (m1 == 2) {
    g1 = 1-s;
  }
  if (m2 == 1) {
    g2 = 1;
  } else if (m2 == 2) {
    g2 = 1-s;
  }
 
  result = g1 / (g1 + g2);

  return result;
}

int main(int argc, char* argv[])
{
  int lattsize, nogen, numruns, numslices;
  double mu, sel, delta_s, delta_mu;

  string statsfilename;  
  ofstream outstats;
  ofstream outstats2;
  ofstream outstats3;
  ifstream testoutstats;
  int filecounter=0;
  string tempstr;

  ostringstream tempstring;	

  MPI::Init(argc, argv);
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  int root = 0;

  int seed = time(0)+rownd(double(rank)*double(time(0))/double(size));

  mt19937 generator(seed);

  if (rank == 0)
  {
    lattsize = atoi(argv[1]);
    nogen = atoi(argv[2]);
    numslices = atoi(argv[3]);
    numruns = atoi(argv[4]);

    tempstring << filecounter;

    statsfilename += "DKphase_run";
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

  }

  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Bcast(&lattsize, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&nogen, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&numruns, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&numslices, 1, MPI::INT, root);

  printf("Rank #%d checking in\n", rank);

  MPI::COMM_WORLD.Barrier();

  mu = 0.02;
  sel = 0.3;

  delta_mu = (0.03 - 0.02)/(double)numslices;

  if (rank == 0) {

    cout << "# Lx: " << lattsize << " Lt: " << nogen << " No. Runs: " << numruns << " No. Slices: " << numslices << endl;
    cout << "{/Symbol m}\ts\tf\t{/Symbol r}_{DW}" << endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " No. Runs: " << numruns << " No. Slices: " << numslices << endl;
    outstats << "{/Symbol m}\ts\tf\t{/Symbol r}_{DW}" << endl;

  }

  vector< int > lattodd(lattsize);
  vector< int > latteven(lattsize);
    
  double dwden, sdwden;

  double f, sf;

  int dwnum;

  for (int m = 0; m <= numslices; m++) {

    f = 0.0;
    sf = 0.0;
    dwden = 0.0;
    sdwden = 0.0;

    for (int q = 1; q <= numruns; q++) {
  
      // initialize lattice
      for (int i = 0; i < lattsize; i++) {
        if (mu > 0.0000000001) {
          latteven[i] = 1;
        } else {
          latteven[i] = rbetween(generator,1,2);
	}
      }

      // population growth
      for (int t=0; t<=nogen; t++) {

        dwnum = 0;

        // odd step
        for (int i = 0; i < lattsize; i++) {
          if (latteven[i] == latteven[mod(i-1,lattsize)]) {
            lattodd[i] = latteven[i];
          } else {
            ++dwnum;
            if (random01(generator) < prob(latteven[i],latteven[mod(i-1,lattsize)],sel)) {
              lattodd[i] = latteven[i];
            } else {
              lattodd[i] = latteven[mod(i-1,lattsize)];
            }
          }
          if (lattodd[i] == 1) {
            if (mu > 0.0000000001) {
              if (random01(generator) < mu) lattodd[i]++;
            }
          }
        }
        ++t;
        if (dwnum == 0 && lattodd[1] == 2) break;
        dwnum = 0;
  
        // even step
        for (int i = 0; i < lattsize; i++) {
          if (lattodd[i] == lattodd[mod(i+1,lattsize)]) {
            latteven[i] = lattodd[i];
          } else {
            ++dwnum;
            if (random01(generator) < prob(lattodd[i],lattodd[mod(i+1,lattsize)],sel)) {
              latteven[i] = lattodd[i];
            } else {
              latteven[i] = lattodd[mod(i+1,lattsize)];
            }
          }
          if (latteven[i] == 1) {
            if (mu > 0.0000000001) {
              if (random01(generator) < mu) latteven[i]++;
            }
          }
        }
        if (dwnum == 0 && latteven[1] == 2) break;

      }
  
      for (int i = 0; i < lattsize; i++) {
       if (lattodd[i] == 1) f += 1.0/(double)lattsize;
      }
  
      dwden += (double)dwnum/(double)lattsize;
    }

    MPI::COMM_WORLD.Barrier();
  
    MPI::COMM_WORLD.Reduce(&f, &sf, 1, MPI::DOUBLE, MPI::SUM, root);

    MPI::COMM_WORLD.Barrier();

    MPI::COMM_WORLD.Reduce(&dwden, &sdwden, 1, MPI::DOUBLE, MPI::SUM, root);

    MPI::COMM_WORLD.Barrier();

    if (rank == 0) {
      outstats << mu << "\t" << sel << "\t" << sf/(double)numruns/(double)size << "\t" << sdwden/(double)numruns/(double)size << endl;
      cout << mu << "\t" << sel << "\t" << sf/(double)numruns/(double)size << "\t" << sdwden/(double)numruns/(double)size << endl;
    }
  mu += delta_mu;
  }

  if (rank == 0) outstats.close();

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize();

  return 0;
  
}
  
  

