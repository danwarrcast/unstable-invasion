#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>
#include <mpi.h>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>

using namespace std;

using boost::mt19937;
using boost::random::uniform_real_distribution;
using boost::random::uniform_int_distribution;

int rownd (double a) {
        return(int(a+ 0.5));
}

bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

// determine winner of cell competition
unsigned int winstrain(unsigned int a, unsigned int b, unsigned int c, unsigned int d, unsigned int e, unsigned int f, vector < vector < vector < vector < vector < vector < vector < double > > > > > > > &gammas, double rannum) {

  double G = 0.0;
  unsigned int result;

  for (int s = 0; s < 6; s++) {
    G += gammas[a][b][c][d][e][f][s];
    if (G >= rannum) {
      if (s == 0) {
        return a;
      } else if (s == 1) {
        return b;
      } else if (s == 2) {
        return c;
      } else if (s == 3) {
        return d;
      } else if (s == 4) {
        return e;
      } else if (s == 5) {
        return f;
      }
    }
  }
  
}

// calls a%b with the result always positive
unsigned int mod(unsigned int a, unsigned int b)
{
  unsigned int r = (a%b+b)%b;
  return r;
}

unsigned int chooseactive(vector< vector < unsigned int > > &alist, vector< unsigned int > &aM, double rN, double sG, double sW) {
  unsigned int index;
  double sum = 0;
  double denom = (1.0-sG)*aM[1] + aM[2] + (1.0-sW)*aM[3];
  vector < double > sums = {0, (1.0-sG)/denom, 1.0/denom, (1.0-sW)/denom};
  for (unsigned int i = 0; i < aM[0]; i++) {
    sum += sums[alist[i][2]];
    if (sum >= rN) {
      index = i;
      break;
    }
  }
  if (sum < rN) cout << "Loop is over and sum is only " << sum << " so I send index " << index << endl;
  return index;
}

int main(int argc, char* argv[])
{
  int lattsize, nogen, numruns;
  double mu, sG, sW;
  double rannum1, rannum2;
  bool image;

  string statsfilename;  
  ofstream outstats, outstats2, outstats3;
  ifstream testoutstats;
  int filecounter=0;
  string tempstr;

  ostringstream tempstring;	

  MPI::Init(argc, argv);
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  int root = 0;

  int seed = time(0)+rownd(double(rank)*double(time(0))/double(size));

  if (rank == 0) {

      lattsize = atoi(argv[4]);
      nogen = atoi(argv[5]);
      numruns = atoi(argv[6]);
      sW = atof(argv[1]);
      sG = sW - atof(argv[2]);
      mu = atof(argv[3]);
      image = false;

    if (argc == 8 && to_bool(argv[7])) { image = true; }
    else if (argc == 8 && !to_bool(argv[7])) { image = false; }

    statsfilename += "diffusion_out/man/ranseq";

    tempstring.str("");
    tempstring.clear();
    tempstring << "_run" << filecounter;
    statsfilename += tempstring.str();

    tempstring.str("");
    tempstring.clear();
    tempstring << filecounter;

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

    statsfilename = "diffusion_out/man/ranseq_w";

    tempstring.str("");
    tempstring.clear();
    tempstring << "_run" << filecounter;
    statsfilename += tempstring.str();

    outstats2.open(statsfilename.c_str());

    cout << statsfilename.c_str() << endl;

    if (image) {
      statsfilename = "diffusion_out/man/ranseq_image";

      tempstring.str("");
      tempstring.clear();
      tempstring << "_run" << filecounter;
      statsfilename += tempstring.str();

      outstats3.open(statsfilename.c_str());

      cout << statsfilename.c_str() << endl;
    }
  }
  
  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Bcast(&filecounter, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&lattsize, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&nogen, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&numruns, 1, MPI::INT, root);
  MPI::COMM_WORLD.Bcast(&sW, 1, MPI::DOUBLE, root);
  MPI::COMM_WORLD.Bcast(&sG, 1, MPI::DOUBLE, root);
  MPI::COMM_WORLD.Bcast(&mu, 1, MPI::DOUBLE, root);

  if (rank == 0) {
    cout << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " N_procs: " << size << " Has_Image: " << image << endl;
    cout << " # sW = " << sW << " b = " << sW - sG << " mu = " << mu <<  endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " N_procs: " << size << " Has_Image: " << image << endl;
    outstats << " # sW = " << sW << " b = " << sW - sG << " mu = " << mu <<  endl;
  }

  int lattsize2 = lattsize*4;
  
  mt19937 generator(seed);
  boost::random::uniform_real_distribution< double > dis(0,1);

    // calculate all gamma factors for competition
  vector< vector< vector< vector< vector< vector< vector< double> > > > > > > gammas(4, vector< vector< vector< vector< vector< vector< double > > > > > >(4, vector< vector< vector< vector< vector< double > > > > >(4, vector< vector< vector< vector< double > > > >(4, vector< vector< vector< double > > >(4, vector< vector< double > >(4, vector< double >(6)))))));
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
        for (int l = 1; l <= 3; l++) {
          for (int m = 1; m <= 3; m++) {
            for (int n = 1; n <= 3; n++) {
              double g_a = 1.0 - sW*(double)(i-1)*(double)(i-2)/2.0 - sG*(double)(2-i)*(double)(3-i)/2.0;
              double g_b = 1.0 - sW*(double)(j-1)*(double)(j-2)/2.0 - sG*(double)(2-j)*(double)(3-j)/2.0;
              double g_c = 1.0 - sW*(double)(k-1)*(double)(k-2)/2.0 - sG*(double)(2-k)*(double)(3-k)/2.0;
              double g_d = 1.0 - sW*(double)(l-1)*(double)(l-2)/2.0 - sG*(double)(2-l)*(double)(3-l)/2.0;
              double g_e = 1.0 - sW*(double)(m-1)*(double)(m-2)/2.0 - sG*(double)(2-m)*(double)(3-m)/2.0;
              double g_f = 1.0 - sW*(double)(n-1)*(double)(n-2)/2.0 - sG*(double)(2-n)*(double)(3-n)/2.0;
              double Sum = g_a + g_b + g_c + g_d + g_e + g_f;
              gammas[i][j][k][l][m][n][0] = g_a/Sum;
              gammas[i][j][k][l][m][n][1] = g_b/Sum;
              gammas[i][j][k][l][m][n][2] = g_c/Sum;
              gammas[i][j][k][l][m][n][3] = g_d/Sum;
              gammas[i][j][k][l][m][n][4] = g_e/Sum;
              gammas[i][j][k][l][m][n][5] = g_f/Sum;
            }
          }
        }
      }
    }
  }

  bool pow_2 = false;
  unsigned int t2;
  unsigned int nogen2 = ceil(log2(nogen));

  int left = (lattsize2)/2 - lattsize/4, right = (lattsize2)/2 + lattsize/4;
  vector< vector< unsigned int > >::iterator Lit;
  vector< vector< unsigned int > >::reverse_iterator Rit;
  unsigned int Lmax, Rmax, LInt, RInt, i_rand, j_rand, active_rand;
  double dt;

  char* _emergencyMemory = new char[16384];
  try { 
    vector< double > results(numruns*(nogen2+1)*lattsize*6);
    vector< vector< vector< unsigned int > > > latt(lattsize, vector< vector< unsigned int > > (lattsize2, vector< unsigned int > (3)));
    vector< vector< unsigned int > > active(lattsize*lattsize2, vector< unsigned int >(3));
    vector< unsigned int > aMeta(4);
  } catch(bad_alloc& ex) {
    delete[] _emergencyMemory;
    cerr << "Not enough memory!" << endl;
    cin.get();
    exit(1);
  }
  vector< double > results(numruns*(nogen2+1)*lattsize*6);
  vector< vector< vector< unsigned int > > > latt(lattsize, vector< vector< unsigned int > > (lattsize2, vector< unsigned int > (3)));
  vector < vector < unsigned int > > active(lattsize*lattsize2, vector < unsigned int >(3));
  vector< unsigned int > aMeta(4);

  for (int q = 0; q < numruns; q++) {  

    unsigned int currindexq = q*(nogen2+1)*lattsize*6;
    aMeta.clear();
    // initialize lattice
    for (unsigned int i = 0; i < lattsize; i++) {
      for (unsigned int j = 0; j < lattsize2; j++) {
        latt[i][j][0] = 1;
        latt[i][j][1] = 0;
        latt[i][j][2] = 0;
        if (j > left && j <= right) {
          latt[i][j][0]++;
        }
        if (j == left || j == left+1 || j == right || j == right+1) {
          active[aMeta[0]] = {i,j,latt[i][j][0]};
          latt[i][j][1] = 1;
          latt[i][j][2] = aMeta[0];
          aMeta[0]++;
          aMeta[latt[i][j][0]]++;
        }
      }
    }

    t2 = 0;
    double t = 0.0;
    while (t <= pow(2,nogen2)+1) {
   
      if (t > 0.0 && log2(t) >= t2) {
        pow_2 = true;
      } else {
        pow_2 = false;
      }

      if (aMeta[0] != aMeta[1]+aMeta[2]+aMeta[3]) cout << "Something wrong with aMeta: (" << aMeta[0] << "," << aMeta[1] << "," << aMeta[2] << "," << aMeta[3] << ")" << endl;
      
      //choose a random lattice site to update from the list of active sites
      rannum1 = dis(generator);
      active_rand = chooseactive(active, aMeta, rannum1, sG, sW);
      i_rand = active[active_rand][0];
      j_rand = active[active_rand][1];
      unsigned int a,b,c,d,e,f,g;
      unsigned int lxindex = mod(i_rand-1,lattsize);
      unsigned int rxindex = mod(i_rand+1,lattsize);
      unsigned int uyindex = mod(j_rand-1,lattsize);
      unsigned int dyindex = mod(j_rand+1,lattsize);
      if ( j_rand % 2 == 0) {
        //find neighbors of the chosen site if j-index is even
        a = latt[lxindex][j_rand][0];
        b = latt[rxindex][j_rand][0];
        c = latt[i_rand][uyindex][0];
        d = latt[rxindex][uyindex][0];
        e = latt[i_rand][dyindex][0];
        f = latt[rxindex][dyindex][0];
        g = latt[i_rand][j_rand][0];
      } else {
        //find neighbors of the chosen site if j-index is odd
        a = latt[lxindex][j_rand][0];
        b = latt[rxindex][j_rand][0];
        c = latt[i_rand][uyindex][0];
        d = latt[lxindex][uyindex][0];
        e = latt[i_rand][dyindex][0];
        f = latt[lxindex][dyindex][0];
        g = latt[i_rand][j_rand][0];
      }
      if ( a == b && b == c && c == d && d == e && e == f ) {
        //if all the neighbors of the chosen site have the same identity, don't pull a random number; replace active site with identity of the neighbors
        if ( a != g ) latt[i_rand][j_rand][0] = a;
      } else {
        //if the neighbors of the chosen cell have different identities, pull a random number and decide the winner to replace the active site
        rannum1 = dis(generator);
        latt[i_rand][j_rand][0] = winstrain(a, b, c, d, e, f, gammas, rannum1);
      }
      if (latt[i_rand][j_rand][0] == 2 && mu > 0.0000000001) {
        //if mutation rate is non-zero, check for a mutation event
        rannum2 = dis(generator);
        if (rannum2 < mu) latt[i_rand][j_rand][0]++;
      }
      active[latt[i_rand][j_rand][2]][2] = latt[i_rand][j_rand][0];
      if (latt[i_rand][j_rand][0] != g) {
        //if the lattice has changed during the last update, update active list
        aMeta[g]--;
        aMeta[latt[i_rand][j_rand][0]]++;
        vector < vector < unsigned int > > tmp(7, vector < unsigned int >(3));
        if (j_rand % 2 == 0) {
          //find neighbors if j-index for current updated cell is even
          tmp = {{lxindex,j_rand,a},
                 {rxindex,j_rand,b},
                 {i_rand,uyindex,c},
                 {rxindex,uyindex,d},
                 {i_rand,dyindex,e},
                 {rxindex,dyindex,f},
                 {i_rand,j_rand,latt[i_rand][j_rand][0]}};
        } else {
          //find neighbors if j-index for current updated cell is odd
          tmp = {{lxindex,j_rand,a},
                 {rxindex,j_rand,b},
                 {i_rand,uyindex,c},
                 {lxindex,uyindex,d},
                 {i_rand,dyindex,e},
                 {lxindex,dyindex,f},
                 {i_rand,j_rand,latt[i_rand][j_rand][0]}};
        }
        for (auto it = tmp.begin(); it != tmp.end(); it++) {
          //for each of the neighbors of the updated cell (including updated cell), check if that cell is now active or inactive
          unsigned int new_i = (*it)[0], new_j = (*it)[1], new_m = (*it)[2];
          bool in_active = latt[new_i][new_j][1];
          unsigned int index = latt[new_i][new_j][2];
          unsigned int ta, tb, tc, td, te, tf;
          unsigned int n_lxindex = mod(new_i-1,lattsize);
          unsigned int n_rxindex = mod(new_i+1,lattsize);
          unsigned int n_uyindex = mod(new_j-1,lattsize);
          unsigned int n_dyindex = mod(new_j+1,lattsize);
          if (new_j % 2 == 0) {
            //find nearest neighbors of site (new_i,new_j) if new_j is even
            ta = latt[n_lxindex][new_j][0];
            tb = latt[n_rxindex][new_j][0];
            tc = latt[new_i][n_uyindex][0];
            td = latt[n_rxindex][n_uyindex][0];
            te = latt[new_i][n_dyindex][0];
            tf = latt[n_rxindex][n_dyindex][0];
          } else { 
            //find nearest neighbors of site (new_i,new_j) if new_j is odd
            ta = latt[n_lxindex][new_j][0];
            tb = latt[n_rxindex][new_j][0];
            tc = latt[new_i][n_uyindex][0];
            td = latt[n_lxindex][n_uyindex][0];
            te = latt[new_i][n_dyindex][0];
            tf = latt[n_lxindex][n_dyindex][0];
          }
          //check if current site is in the active list
          if (in_active) {
            //if current site is in the active list, check if it *should* be in the active list
            if ( new_m == ta && new_m == tb && new_m == tc && new_m == td && new_m == te && new_m == tf) {
              //if current site should not be in active list anymore (has same identity as all its neihbors), remove it
              active[index] = active[aMeta[0]-1];
              latt[active[index][0]][active[index][1]][2] = index;
              latt[new_i][new_j][1] = 0;
              latt[new_i][new_j][2] = active.size()+1;
              aMeta[0]--;
              aMeta[new_m]--;
            }
          } else {
            if (new_m != ta || new_m != tb || new_m != tc || new_m != td || new_m != te || new_m != tf) {
              //if current lattice site is not in the active list and has at least one neighbor with a different identity, add it to active list
              active[aMeta[0]] = {new_i, new_j, new_m};
              latt[new_i][new_j][1] = 1;
              latt[new_i][new_j][2] = aMeta[0];
              aMeta[0]++;
              aMeta[new_m]++;
            }
          }
        }
      }

      if (pow_2) {

        unsigned int currindexqt2 = currindexq + t2*lattsize*6;
        if (rank == 0) cout << rank << ": (q,t2,t) = (" << q << "," << t2 << "," << t << "); aSize = " << aMeta[0] << endl;

        for (int i = 0; i < lattsize; i++) {

          //Find first and last elements in current column of the lattice with value >1
          for (auto it = latt[i].begin(); it != latt[i].end(); it++) {
            if ((*it)[0] > 1) {
              Lit = it;
              break;
            }
          }
          for (auto it = latt[i].rbegin(); it != latt[i].rend(); it++) {
            if ((*it)[0] > 1) {
              Rit = it;
              break;
            }
          }
          //The first and last elements in the current column with value >1 will act as our maximum range for this column
          //since the the identity for j<Lmax and j>Ramax must be 1
          Lmax = distance(latt[i].begin(),Lit);
          Rmax = distance(latt[i].begin(),(Rit+1).base());

          std::vector< int > Lmatches;
          std::vector< int > Rmatches;
          //Find all green cells (M = 1) within the range Lmax < j < Rmax and add these cells to the vector Lmatches or Rmatches
          //depending on which side of the lattice the green cell is on
          for (auto j = Lit, toofar = (Rit+1).base()+1; j != toofar; ++j) {
            if ( (*j)[0] == 1 && distance(latt[i].begin(),j) < lattsize2/2 ) Lmatches.push_back(distance(latt[i].begin(),j));
              if ( (*j)[0] == 1 && distance(latt[i].begin(),j) > lattsize2/2 ) Rmatches.push_back(distance(latt[i].begin(),j));
          }

          double Lavg = 0, Ravg = 0;
          //average the j-indices of the green cells for each column to determine the average "position" of the boundary in this column
          if (Lmatches.size() > 0) {
            LInt = Lmatches[Lmatches.size() - 1] - Lmax;
            Lavg = (Lmax - 1.0)/( (double)Lmatches.size() + 1.0 );
            for (auto k = Lmatches.begin(); k != Lmatches.end(); ++k) {
              Lavg += *k/( (double)Lmatches.size() + 1.0 );
            }
            if (Lavg < 0.000001) Lavg = Lmax;
          } else {
            Lavg = Lmax;
            LInt = 0;
          }
          if (Rmatches.size() > 0) {
            RInt = Rmax - Rmatches[0];
            Ravg = (Rmax + 1.0)/( (double)Rmatches.size() + 1.0 );
            for (auto k = Rmatches.begin(); k != Rmatches.end(); ++k) {
              Ravg += *k/( (double)Rmatches.size() + 1.0 );
            }
            if (Ravg < 0.000001) Ravg = Rmax;
          } else {
            Ravg = Rmax;
            RInt = 0;
          }
          int currindex = currindexqt2+i*6;
          results[currindex+0] = Lavg;
          results[currindex+1] = Ravg;
          results[currindex+2] = Lmax;
          results[currindex+3] = Rmax;
          results[currindex+4] = LInt;
          results[currindex+5] = RInt;
        }
        pow_2 = false;
        t2++;
      }
      //Advance the time according to Gillespie's algorithm
      dt = -log(1.0 - dis(generator))/(double)aMeta[0];
      t = t + dt;
    }
  }

  if (rank == 0 && image) {
    for (auto it = latt.begin(); it != latt.end(); ++it) {
      for (auto jt = (*it).begin(); jt != (*it).end(); ++jt) {
        double x = (double)distance(latt.begin(), it)-0.5*mod(distance((*it).begin(), jt),2);
        double y = (double)distance((*it).begin(), jt)*sqrt(3)/2;
        if ((*jt)[0] > 1) outstats3 << t2 << "\t" << x << "\t" << y << "\t" <<  (*jt)[0] << "\t" << (*jt)[1] << endl;
      }
    }
    outstats3.close();
  }

  //delete latt and active vectors to clear memory for new vectors for analysis
  latt.clear();
  vector< vector< vector< unsigned int > > >().swap(latt);
  active.clear();
  vector< vector< unsigned int > >().swap(active);

  vector< vector< vector< double > > > avgs(numruns, vector< vector< double > > (nogen2+1, vector< double > (6)));
  vector< vector< vector< double > > > avgs2(numruns, vector< vector< double > > (nogen2+1, vector< double > (6)));

  for (int i = 0; i < numruns; i++) {
    int currindexq = i*(nogen2+1)*lattsize*6;
    for (int j = 0; j <= nogen2; j++) {
      int currindexqt2 = currindexq+j*lattsize*6;
      for (int k = 0; k < lattsize; k++) {
        int currindex = currindexqt2 + k*6;
        for (int l = 0; l < 6; l++) avgs[i][j][l] += results[currindex+l]/(double)lattsize;
      }
    }
  }
  for (int i = 0; i < numruns; i++) {
    int currindexq = i*(nogen2+1)*lattsize*6;
    for (int j = 0; j <= nogen2; j++) {
      int currindexqt2 = currindexq + j*lattsize*6;
      for (int k = 0; k < lattsize; k++) {
        int currindex = currindexqt2 + k*6;
        for (int l = 0; l < 6; l++) avgs2[i][j][l] += pow(avgs[i][j][l] - results[currindex+l],2)/(double)lattsize;
      }
    }
  }

  //clear results vector to clear memory for new vectors
  results.clear();
  vector< double >().swap(results);
 
  vector< double > w_avgs((nogen2+1)*3);
  vector< double > sw_avgs((nogen2+1)*3);

  for (int i = 0; i < numruns; i++) {
    for (int j = 0; j <= nogen2; j++) {
      for (int k = 0; k < 3; k++) w_avgs[j*3+k] += (avgs2[i][j][2*k] + avgs2[i][j][2*k+1])/2.0/(double)numruns;
    }
  }

  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Reduce(&w_avgs.front(), &sw_avgs.front(), (nogen2+1)*3, MPI::DOUBLE, MPI::SUM, root);

  MPI::COMM_WORLD.Barrier();

  if (rank == 0) {
    
    //only print results vector for run q = 0, t2 = nogen2
    int currindext2 = nogen2*lattsize*6;
    for (int k = 0; k < lattsize; k++) {
      int currindex = currindext2+k*6;
      outstats << "0\t" << nogen2 << "\t" << k << "\t" << results[currindex+0]/(double)size
                                            << "\t" << results[currindex+1]/(double)size
                                            << "\t" << results[currindex+2]/(double)size
                                            << "\t" << results[currindex+3]/(double)size
                                            << "\t" << results[currindex+4]/(double)size
                                            << "\t" << results[currindex+5]/(double)size << endl;
    }
    for (int i = 0; i <= nogen2; i++) outstats2 << pow(2,i) << "\t" << sqrt(sw_avgs[i*3+0]/(double)size)
                                                           << "\t" << sqrt(sw_avgs[i*3+1]/(double)size)
                                                           << "\t" << sqrt(sw_avgs[i*3+2]/(double)size) << endl;
  }

  MPI::COMM_WORLD.Barrier();

  outstats.close();

  outstats2.close();

  MPI::COMM_WORLD.Barrier();

  MPI::Finalize();

  return 0;
  
}
