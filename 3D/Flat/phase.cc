#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>
#include <mpi.h>

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
  int lattsize, nogen, numruns, numslices;
  double mu, sG, sW, delta_s, delta_mu;
  double rannum1, rannum2;

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

  int seed = time(0);

  mt19937 generator(seed);
  boost::random::uniform_real_distribution< double > dis(0,1);
  boost::random::uniform_int_distribution< int > disi12(1,2);
  boost::random::uniform_int_distribution< int > disi13(1,3);

  if (rank == 0) {
    lattsize = atoi(argv[1]);
    nogen = atoi(argv[2]);
    numslices = atoi(argv[3]);
    numruns = atoi(argv[4]);

    tempstring << filecounter;

    statsfilename += "3Dphase_run";
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

  sG = 0.0;
  sW = 0.3;
  mu = 0.0;
  delta_s = (sW)/(double)numslices;
  delta_mu = (0.04)/(double)numslices;

  if (rank == 0) {
    cout << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << size*numruns << " N_slices: " << numslices << endl;
    cout << " # sW = " << sW << endl;
    cout << "{/Symbol m}\ts_{G}\tf_{G}\tf_{W}" << endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << size*numruns << " N_slices: " << numslices << endl;
    outstats << " # sW = " << sW << endl;
    outstats << "{/Symbol m}\ts_{G}\tf_{G}\tf_{W}" << endl;
  }

  vector< vector< int > > lattup(lattsize, vector< int >(lattsize));
  vector< vector< int > > lattdown(lattsize, vector< int >(lattsize));

  double fG, fW, sfG, sfW;

  double mu_i = mu, sG_i = sG;

  int num_B, num_M;

  while (mu <= ((numslices + 1) * delta_mu)) {

    if (mu > mu_i) sG = 0.0;

    while (sG <= ((numslices + 1) * delta_s)) {

      fG = 0.0;
      fW = 0.0;
      sfG = 0.0;
      sfW = 0.0;

      for (int q = 1; q <= numruns; q++) {  
  
        // initialize lattice
        for (int i = 0; i < lattsize; i++) {
          for (int j = 0; j < lattsize; j++) {
            if (mu > 0.0000000001) {
              lattdown[i][j] = disi12(generator);
            } else {
              lattdown[i][j] = disi13(generator);
            }
          }
        }

        for (int t=1; t<=nogen; t++) {

          num_B = 0;
          num_M = 0;
          
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
              if (lattup[i][j] == 1) {
                num_B++;
              } else {
                num_M++;
              }
            }
          }
          // end odd step
          if (num_B == 0 || num_M == 0) break;
          ++t;

          num_B = 0;
          num_M = 0; 
 
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
              if (lattup[i][j] == 1) {
                num_B++;
              } else {
                num_M++;
              }
            }
          }
          // end even step
          if (num_B == 0 || num_M == 0) break;
        }
        for (int i = 0; i < lattsize; i++) {
          for (int j = 0; j < lattsize; j++) {
            if (lattdown[i][j] == 1) fG += 1.0/pow((double)lattsize,2);
            if (lattdown[i][j] == 3) fW += 1.0/pow((double)lattsize,2);
          }
        }
      }

      MPI::COMM_WORLD.Barrier();
 
      MPI::COMM_WORLD.Reduce(&fG, &sfG, 1, MPI::DOUBLE, MPI::SUM, root);

      MPI::COMM_WORLD.Barrier();

      MPI::COMM_WORLD.Reduce(&fW, &sfW, 1, MPI::DOUBLE, MPI::SUM, root);

      MPI::COMM_WORLD.Barrier();

      if (rank == 0) {
        outstats << mu << "\t" << sG << "\t" << sfG/(double)numruns/(double)size << "\t" << sfW/(double)numruns/(double)size << endl;
        cout << mu << "\t" << sG << "\t" << sfG/(double)numruns/(double)size << "\t" << sfW/(double)numruns/(double)size << endl;
      }
      sG += delta_s;
    }
    mu += delta_mu;
  }

  if (rank == 0) outstats.close();

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize();

  return 0;
  
}
  
  

