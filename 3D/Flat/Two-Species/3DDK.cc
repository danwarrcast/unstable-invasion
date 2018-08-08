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
int winstrain(int a, int b, int c, double sel, double rannum) {

  double g_a = 1.0 - sel*(double)(a-1);
  double g_b = 1.0 - sel*(double)(b-1);
  double g_c = 1.0 - sel*(double)(c-1);

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
  double mu, sel, delta_s, delta_mu;
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

  if (rank == 0) {
    lattsize = atoi(argv[1]);
    nogen = atoi(argv[2]);
    numslices = atoi(argv[3]);
    numruns = atoi(argv[4]);

    tempstring << filecounter;

    statsfilename += "3DDK_run";
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

  sel = 0.0;
  mu = 0.0;
  delta_s = (0.4)/(double)numslices;
  delta_mu = (0.2)/(double)numslices;

  if (rank == 0) {
    cout << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << size*numruns << " N_slices: " << numslices << endl;
    cout << "{/Symbol m}\ts_{G}\tf_{R}\tf_{B}" << endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << size*numruns << " N_slices: " << numslices << endl;
    outstats << "{/Symbol m}\ts_{G}\tf_{R}\tf_{B}" << endl;
  }

  vector< vector< int > > lattup(lattsize, vector< int >(lattsize));
  vector< vector< int > > lattdown(lattsize, vector< int >(lattsize));

  double fR, fB, sfR, sfB;

  double mu_i = mu, sel_i = sel;

  int num_R, num_B;

  while (mu <= ((numslices + 1) * delta_mu)) {

    if (mu > mu_i) sel = 0.0;

    while (sel <= ((numslices + 1) * delta_s)) {

      fR = 0.0;
      fB = 0.0;
      sfR = 0.0;
      sfB = 0.0;

      for (int q = 1; q <= numruns; q++) {  
  
        // initialize lattice
        for (int i = 0; i < lattsize; i++) {
          for (int j = 0; j < lattsize; j++) {
            if (mu > 0.0000000001) {
              lattdown[i][j] = 1;
            } else {
              lattdown[i][j] = disi12(generator);
            }
          }
        }

        for (int t=1; t<=nogen; t++) {

          num_R = 0;
          num_B = 0;
          
          // odd step
          for (int i = 0; i < lattsize; i++) {
            for (int j = 0; j < lattsize; j++) {

              if ( j % 2 == 0) {
                if ((lattdown[i][j] == lattdown[i][mod(j+1,lattsize)]) && (lattdown[i][j] == lattdown[mod(i+1,lattsize)][mod(j+1,lattsize)])) {
                  lattup[i][j] = lattdown[i][j];
                } else {
                  rannum1 = dis(generator);
                  lattup[i][j] = winstrain(lattdown[i][j], lattdown[i][mod(j+1,lattsize)], lattdown[mod(i+1,lattsize)][mod(j+1,lattsize)], sel, rannum1);
                }
              } else {
                if ((lattdown[i][j] == lattdown[mod(i-1,lattsize)][mod(j+1,lattsize)]) && (lattdown[i][j] == lattdown[i][mod(j+1,lattsize)])) {
                  lattup[i][j] = lattdown[i][j];
                } else {
                  rannum1 = dis(generator);
                  lattup[i][j] = winstrain(lattdown[i][j], lattdown[mod(i-1,lattsize)][mod(j+1,lattsize)], lattdown[i][mod(j+1,lattsize)], sel, rannum1);
                }
              }
              if (lattup[i][j] == 1 && mu > 0.0000000001) {
                rannum2 = dis(generator);
                if (rannum2 < mu) lattup[i][j]++;
              }
              if (lattup[i][j] == 1) {
                num_R++;
              } else {
                num_B++;
              }
            }
          }
          // end odd step
          if (num_R == 0) break;
          ++t;

          num_R = 0;
          num_B = 0; 
 
          // even step
          for (int i = 0; i < lattsize; i++) {
            for (int j = 0; j < lattsize; j++) {

              if ( j % 2 == 0) {
                if ((lattup[i][j] == lattup[i][mod(j-1,lattsize)]) && (lattup[i][j] == lattup[mod(i+1,lattsize)][mod(j-1,lattsize)])) {
                  lattdown[i][j] = lattup[i][j];
                } else {
                  rannum1 = dis(generator);
                  lattdown[i][j] = winstrain(lattup[i][j], lattup[i][mod(j-1,lattsize)], lattup[mod(i+1,lattsize)][mod(j-1,lattsize)], sel, rannum1);
                }
              } else {
                if ((lattup[i][j] == lattup[i][mod(j-1,lattsize)]) && (lattup[i][j] == lattup[mod(i-1,lattsize)][mod(j-1,lattsize)])) {
                  lattdown[i][j] = lattup[i][j];
                } else {
                  rannum1 = dis(generator);
                  lattdown[i][j] = winstrain(lattup[i][j], lattup[i][mod(j-1,lattsize)], lattup[mod(i-1,lattsize)][mod(j-1,lattsize)], sel, rannum1);
                }
              }
              if (lattdown[i][j] == 1 && mu > 0.0000000001) {
                rannum2 = dis(generator);
                if (rannum2 < mu) lattdown[i][j]++;
              }
              if (lattup[i][j] == 1) {
                num_R++;
              } else {
                num_B++;
              }
            }
          }
          // end even step
          if (num_R == 0) break;
        }
        for (int i = 0; i < lattsize; i++) {
          for (int j = 0; j < lattsize; j++) {
            if (lattdown[i][j] == 1) fR += 1.0/pow((double)lattsize,2);
            if (lattdown[i][j] == 2) fB += 1.0/pow((double)lattsize,2);
          }
        }
      }

      MPI::COMM_WORLD.Barrier();
 
      MPI::COMM_WORLD.Reduce(&fR, &sfR, 1, MPI::DOUBLE, MPI::SUM, root);

      MPI::COMM_WORLD.Barrier();

      MPI::COMM_WORLD.Reduce(&fB, &sfB, 1, MPI::DOUBLE, MPI::SUM, root);

      MPI::COMM_WORLD.Barrier();

      if (rank == 0) {
        outstats << mu << "\t" << sel << "\t" << sfR/(double)numruns/(double)size << "\t" << sfB/(double)numruns/(double)size << endl;
        cout << mu << "\t" << sel << "\t" << sfR/(double)numruns/(double)size << "\t" << sfB/(double)numruns/(double)size << endl;
      }
      sel += delta_s;
    }
    mu += delta_mu;
  }

  if (rank == 0) outstats.close();

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize();

  return 0;
  
}
  
  

