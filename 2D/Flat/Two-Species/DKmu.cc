#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>
#include <mpi.h>

using namespace std;

using boost::mt19937;

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
  
  g1 = 1.0 - s*(m1-1);
  g2 = 1.0 - s*(m2-1);

  result = g1 / (g1 + g2);

  return result;
}

int main(int argc, char* argv[])
{
  int lattsize, nogen, numruns, numslices;
  double mu, sel, delta_mu;
  double rannum1, rannum2;

  string statsfilename;  
  ofstream outstats;
  ofstream outstats_f;
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
  boost::random::uniform_real_distribution< double > dis(0,1);
  boost::random::uniform_int_distribution< int > disi(1,2);

  if (rank == 0)
  {
    lattsize = atoi(argv[1]);
    nogen = atoi(argv[2]);
    numslices = atoi(argv[3]);
    numruns = atoi(argv[4]);

    tempstring << filecounter;

    statsfilename += "DKmu_run";
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

    statsfilename.erase(statsfilename.end()-tempstr.size()-4,statsfilename.end());
    statsfilename += "f_run" + tempstring.str();
    
    outstats_f.open(statsfilename.c_str());

    cout << statsfilename.c_str() << endl;

  }

  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Bcast(&lattsize, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&nogen, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&numruns, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&numslices, 1, MPI::INT, root);

  printf("Rank #%d checking in\n", rank);

  MPI::COMM_WORLD.Barrier();

  mu = 0.0275;
  sel = 0.3;

  delta_mu = (0.028-0.0275)/(double)numslices;

  if (rank == 0) {

    cout << "# Lx: " << lattsize << " Lt: " << nogen << " No. Runs: " << numruns << " No. Slices: " << numslices << endl;
    cout << "# s: " << sel << endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " No. Runs: " << numruns << " No. Slices: " << numslices << endl;
    outstats << "# s: " << sel << endl;

    outstats_f << "# Lx: " << lattsize << " Lt: " << nogen << " No. Runs: " << numruns << " No. Slices: " << numslices << endl;
    outstats_f << "# s: " << sel << endl;

  }

  vector< int > lattodd(lattsize);
  vector< int > latteven(lattsize);


  vector< vector< double > > dwden(numslices+1, vector< double >(floor(log2(nogen))+1, 0.0));
  vector< vector< double > > sdwden(numslices+1, vector< double >(floor(log2(nogen))+1, 0.0));

  vector< vector< double >  > f(numslices+1, vector< double >(floor(log2(nogen))+1, 0.0));
  vector< vector< double > > sf(numslices+1, vector< double >(floor(log2(nogen))+1, 0.0));

  int count=0, t2;
  bool is_pow_2;
  double mu_i = mu;

  for (int m = 0; m <= numslices; m++) {

    for (int q = 1; q <= numruns; q++) {
  
      // initialize lattice
      for (int i = 0; i < lattsize; i++) {
        if (mu > 0.0000000001) {
          latteven[i] = 1;
        } else {
          latteven[i] = disi(generator);
	}
      }

      count = 1;
      t2 = 0;

      // population growth
      for (int t=0; t<=nogen; t++) {

	is_pow_2 = false;

	if ( t == count ) {
	  is_pow_2 = true;
	  if (t > 1) t2++;
	  count = count*2;
	}

        // odd step
        for (int i = 0; i < lattsize; i++) {
          if (latteven[i] == latteven[mod(i-1,lattsize)]) {
            lattodd[i] = latteven[i];
          } else {
            if (is_pow_2) dwden[m][t2] += 1.0/(double)lattsize;

	    rannum1 = dis(generator);
            if (rannum1 < prob(latteven[i],latteven[mod(i-1,lattsize)],sel)) {
              lattodd[i] = latteven[i];
            } else {
              lattodd[i] = latteven[mod(i-1,lattsize)];
            }
          }
          if (lattodd[i] == 1) {
            if (mu > 0.0000000001) {
	      rannum2 = dis(generator);
              if (rannum2 < mu) lattodd[i]++;
            }
          }
	  if (is_pow_2) {
            if (lattodd[i] == 1) f[m][t2] += 1.0/(double)lattsize;
	  }
	}
        ++t;
	is_pow_2 = false;

	if ( t == count ) {
	  is_pow_2 = true;
	  if (t > 1) t2++;
	  count = count*2;
	}

        // even step
        for (int i = 0; i < lattsize; i++) {
          if (lattodd[i] == lattodd[mod(i+1,lattsize)]) {
            latteven[i] = lattodd[i];
          } else {
            if (is_pow_2) dwden[m][t2] += 1.0/(double)lattsize;

	    rannum1 = dis(generator);
            if (rannum1 < prob(lattodd[i],lattodd[mod(i+1,lattsize)],sel)) {
              latteven[i] = lattodd[i];
            } else {
              latteven[i] = lattodd[mod(i+1,lattsize)];
            }
          }
          if (latteven[i] == 1) {
            if (mu > 0.0000000001) {
	      rannum2 = dis(generator);
              if (rannum2 < mu) latteven[i]++;
            }
          }
	  if (is_pow_2) {
            if (latteven[i] == 1) f[m][t2] += 1.0/(double)lattsize;
	  }
        }

      }
  
    }
  mu += delta_mu;
  }

  MPI::COMM_WORLD.Barrier();

  for (int i = 0; i <= numslices; i++) {
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Reduce(&dwden[i].front(), &sdwden[i].front(), floor(log2(nogen))+1, MPI::DOUBLE, MPI::SUM, root);
    MPI::COMM_WORLD.Reduce(&f[i].front(), &sf[i].front(), floor(log2(nogen))+1, MPI::DOUBLE, MPI::SUM, root);
    MPI::COMM_WORLD.Barrier();
  }

  MPI::COMM_WORLD.Barrier();

  if (rank == 0) {
    outstats << "mu:";
    outstats_f << "mu:";
    for (int j = 0; j <= numslices; j++) outstats << "\t" << mu_i + j*delta_mu; 
    for (int j = 0; j <= numslices; j++) outstats_f << "\t" << mu_i + j*delta_mu; 
    outstats << endl;
    outstats_f << endl;
    for (int i = 0; i < floor(log2(nogen))+1; i++) {
      outstats << pow(2,i);
      outstats_f << pow(2,i);
      for (int j = 0; j <= numslices; j++) outstats << "\t" << sdwden[j][i]/(double)numruns/(double)size;
      for (int j = 0; j <= numslices; j++) outstats_f << "\t" << sf[j][i]/(double)numruns/(double)size;
      outstats << endl;
      outstats_f << endl;
    }
  }
 
  MPI::COMM_WORLD.Barrier();

  if (rank == 0) outstats.close();

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize();

  return 0;
  
}
  
  

