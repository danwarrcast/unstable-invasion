#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>
#include <mpi.h>

using namespace std;

using boost::mt19937;
using boost::random::uniform_real_distribution;

int rownd (double a) {
        return(int(a+ 0.5));
}

// calls a%b with the result always positive
int mod(int a, int b)
{
  int r = a%b;
  return r<0? r+b:r;
}

// determine probability of strain m1 winning over strain m2
double prob(int m1, int m2, double DPs, double Bys)
{
  double g1, g2, result;

  g1 = 1.0-Bys*(double)((3-m1)*(2-m1))/2.0-DPs*(double)((m1-2)*(m1-1))/2.0;
  g2 = 1.0-Bys*(double)((3-m2)*(2-m2))/2.0-DPs*(double)((m2-2)*(m2-1))/2.0;

  result = g1 / (g1 + g2);

  return result;
}

int main(int argc, char* argv[])
{
  int lattsize, nogen, numruns;
  double mu, DPsel, Bysel;
  double rannum1, rannum2;

  string statsfilename;  
  ofstream outstats;
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
  uniform_real_distribution< double > dis(0,1);
  uniform_int_distribution< int > disi(1,2);

  if (rank == 0)
  {
    lattsize = atoi(argv[1]);
    nogen = atoi(argv[2]);
    numruns = atoi(argv[3]);

    tempstring << filecounter;

    statsfilename += "diffusion_run";
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

  MPI::COMM_WORLD.Barrier();

  DPsel = 0.3;
  Bysel = 0.245;
  mu = 0.0276;


  if (rank == 0) {

    cout << "# Lx: " << lattsize << " Lt: " << nogen << " Nruns: " << size*numruns << endl;
    cout << "# sW: " << DPsel << " sG: " << Bysel << " mu: " << mu << endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " Nruns: " << size*numruns << endl;
    outstats << "# sW: " << DPsel << " sG: " << Bysel << " mu: " << mu << endl;
    outstats << "t\tD(t)" << endl;

  }

  vector< int > lattodd(lattsize);
  vector< int > latteven(lattsize);

  // round up total time to the next power of 2  
  nogen = pow(2,ceil(log2(nogen)));

  vector< double > DW(ceil(log2(nogen))+1);
  vector< double > sDW(ceil(log2(nogen))+1);
  vector< double > DW2(ceil(log2(nogen))+1);
  vector< double > sDW2(ceil(log2(nogen))+1);

  bool pow_2;
  int count, t2;
  
  int left = lattsize/4, right = 3*lattsize/4;

  for (int q = 0; q < numruns; q++) {

    // initialize lattice
    for (int i = 0; i < lattsize; i++) {
      latteven[i] = 1;
      if (i > left && i <= right) {
        latteven[i]++;
      }
    }

    count = 1;
    t2 = 0;

    // population growth
    for (int t=1; t<=nogen; t++) {

      if ( t == count ) {
        pow_2 = true;
        if ( t > 1 ) t2++;
	count = count*2;
      } else {
        pow_2 = false;
      }

      // odd step
      for (int i = 0; i < lattsize; i++) {
	
        if (latteven[i] == latteven[mod(i-1,lattsize)]) {

          lattodd[i] = latteven[i];

        } else {

          if (pow_2 && latteven[mod(i-1,lattsize)] == 1) {
            DW[t2] += (i-1-left)/2.0/(double)numruns;
            DW2[t2] += pow(i-1-left,2)/2.0/(double)numruns;
	  }
	  if (pow_2 && latteven[i] == 1) {
            DW[t2] += (i-1-right)/2.0/(double)numruns;
            DW2[t2] += pow(i-1-right,2)/2.0/(double)numruns;
	  }

	  rannum1 = dis(generator);
          if (rannum1 < prob(latteven[i],latteven[mod(i-1,lattsize)],DPsel,Bysel)) {
            lattodd[i] = latteven[i];
          } else {
            lattodd[i] = latteven[mod(i-1,lattsize)];
          }

        }
        if (lattodd[i] == 2) {
          rannum2 = dis(generator);
          if (rannum2 < mu) lattodd[i]++; 
        }
      }

      ++t;

      if ( t == count ) {
        pow_2 = true;
        if ( t > 1 ) t2++;
	count = count*2;
      } else {
        pow_2 = false;
      }

      // even step
      for (int i = 0; i < lattsize; i++) {

        if (lattodd[i] == lattodd[mod(i+1,lattsize)]) {

          latteven[i] = lattodd[i];

        } else {

          if (pow_2 && lattodd[i] == 1) {
            DW[t2] += (i-left)/2.0/(double)numruns;
            DW2[t2] += pow(i-left,2)/2.0/(double)numruns;
	  }
	  if (pow_2 && lattodd[mod(i+1,lattsize)] == 1) {
            DW[t2] += (i-right)/2.0/(double)numruns;
            DW2[t2] += pow(i-right,2)/2.0/(double)numruns;
	  }

          rannum1 = dis(generator);
          if (rannum1 < prob(lattodd[i],lattodd[mod(i+1,lattsize)],DPsel,Bysel)) {
            latteven[i] = lattodd[i];
          } else {
            latteven[i] = lattodd[mod(i+1,lattsize)];
          }

        }

        if (latteven[i] == 2) {
	  rannum2 = dis(generator);
          if (rannum2 < mu) latteven[i]++;
	}
      }
    }
  }

  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Reduce(&DW.front(), &sDW.front(), ceil(log2(nogen))+1, MPI::DOUBLE, MPI::SUM, root);

  MPI::COMM_WORLD.Reduce(&DW2.front(), &sDW2.front(), ceil(log2(nogen))+1, MPI::DOUBLE, MPI::SUM, root);

  MPI::COMM_WORLD.Barrier();

  if (rank == 0) {
    for (int i = 0; i < sDW.size(); i++) {
      outstats << pow(2,i) << "\t" << sqrt(sDW2[i]/(double)size - pow(sDW[i]/(double)size,2)) << endl;
    }
    outstats.close();
  }

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize();

  return 0;
  
}
  
  

