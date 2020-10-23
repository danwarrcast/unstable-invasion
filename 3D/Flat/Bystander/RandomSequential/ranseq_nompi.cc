#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <math.h>
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
using boost::random::uniform_int_distribution;
using boost::random::uniform_real_distribution;

int rownd(double a)
{
  return (int(a + 0.5));
}

bool to_bool(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

// determine winner of cell competition
int winstrain(int a, int b, int c, int d, int e, int f, vector<vector<vector<vector<vector<vector<vector<double>>>>>>> &gammas, double rannum)
{

  double G = 0.0;
  int winner;

  for (int s = 0; s < 4; s++)
  {
    G += gammas[a][b][c][d][e][f][s];
    if (G >= rannum)
    {
      winner = s;
      break;
    }
  }
  return winner;
}

// calls a%b with the result always positive
int mod(int a, int b)
{
  int r = (a % b + b) % b;
  return r;
}

int chooseactive(vector<vector<int>> &alist, vector<int> &aM, double rN, double sG, double sW)
{
  int index;
  double sum = 0;
  double denom = (1.0 - sG) * aM[1] + aM[2] + (1.0 - sW) * aM[3];
  vector<double> sums = {0, (1.0 - sG) / denom, 1.0 / denom, (1.0 - sW) / denom};
  for (int i = 0; i < aM[0]; i++)
  {
    sum += sums[alist[i][2]];
    if (sum >= rN)
    {
      index = i;
      break;
    }
  }
  if (sum < rN)
    cout << "Loop is over and sum is only " << sum << " so I send index " << index << endl;
  return index;
}

int main(int argc, char *argv[])
{
  int lattsize, nogen, numruns;
  double mu, sG, sW;
  double rannum1, rannum2;
  bool image;

  string statsfilename;
  ofstream outstats, outstats2, outstats3;
  ifstream testoutstats;
  int filecounter = 0;
  string tempstr;

  ostringstream tempstring;

  int seed = time(0);

  lattsize = atoi(argv[4]);
  nogen = atoi(argv[5]);
  numruns = atoi(argv[6]);
  sW = atof(argv[1]);
  sG = sW - atof(argv[2]);
  mu = atof(argv[3]);
  image = false;

  if (argc == 8 && to_bool(argv[7]))
  {
    image = true;
  }
  else if (argc == 8 && !to_bool(argv[7]))
  {
    image = false;
  }

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
    statsfilename.erase(statsfilename.end() - tempstr.size(), statsfilename.end());
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

  if (image)
  {
    statsfilename = "diffusion_out/man/ranseq_image";

    tempstring.str("");
    tempstring.clear();
    tempstring << "_run" << filecounter;
    statsfilename += tempstring.str();

    outstats3.open(statsfilename.c_str());

    cout << statsfilename.c_str() << endl;
  }

  cout << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " Has_Image: " << image << endl;
  cout << " # sW = " << sW << " b = " << sW - sG << " mu = " << mu << endl;

  outstats << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " Has_Image: " << image << endl;
  outstats << " # sW = " << sW << " b = " << sW - sG << " mu = " << mu << endl;

  int lattsize2 = lattsize * 4;

  mt19937 generator(seed);
  boost::random::uniform_real_distribution<double> dis(0, 1);

  // calculate all gamma factors for competition
  vector<vector<vector<vector<vector<vector<vector<double>>>>>>> gammas(4, vector<vector<vector<vector<vector<vector<double>>>>>>(4, vector<vector<vector<vector<vector<double>>>>>(4, vector<vector<vector<vector<double>>>>(4, vector<vector<vector<double>>>(4, vector<vector<double>>(4, vector<double>(4)))))));
  for (int i = 1; i <= 3; i++)
  {
    for (int j = 1; j <= 3; j++)
    {
      for (int k = 1; k <= 3; k++)
      {
        for (int l = 1; l <= 3; l++)
        {
          for (int m = 1; m <= 3; m++)
          {
            for (int n = 1; n <= 3; n++)
            {
              vector<int> tmp = {0, 0, 0, 0};
              tmp[i]++;
              tmp[j]++;
              tmp[k]++;
              tmp[l]++;
              tmp[m]++;
              tmp[n]++;
              double Denom = (1.0 - sG) * (double)tmp[1] + (double)tmp[2] + (1.0 - sW) * (double)tmp[3];
              gammas[i][j][k][l][m][n][1] = (1.0 - sG) * (double)tmp[1] / Denom;
              gammas[i][j][k][l][m][n][2] = (double)tmp[2] / Denom;
              gammas[i][j][k][l][m][n][3] = (1.0 - sW) * (double)tmp[3] / Denom;
            }
          }
        }
      }
    }
  }

  bool pow_2 = false;
  int t2;
  int nogen2 = ceil(log2(nogen));

  int left = (lattsize2) / 2 - lattsize / 4, right = (lattsize2) / 2 + lattsize / 4;
  vector<vector<int>>::iterator Lit;
  vector<vector<int>>::reverse_iterator Rit;
  int Lmax, Rmax, LInt, RInt, i_rand, j_rand, active_rand;
  double dt;

  char *_emergencyMemory = new char[16384];
  try
  {
    vector<double> results(numruns * (nogen2 + 1) * lattsize * 6);
    vector<vector<vector<int>>> latt(lattsize, vector<vector<int>>(lattsize2, vector<int>(3)));
    vector<vector<int>> active(lattsize * lattsize2, vector<int>(3));
    vector<int> aMeta(4);
  }
  catch (bad_alloc &ex)
  {
    delete[] _emergencyMemory;
    cerr << "Not enough memory!" << endl;
    cin.get();
    exit(1);
  }
  vector<double> results(numruns * (nogen2 + 1) * lattsize * 6);
  vector<vector<vector<int>>> latt(lattsize, vector<vector<int>>(lattsize2, vector<int>(3)));
  vector<vector<int>> active(lattsize * lattsize2, vector<int>(3));
  vector<int> aMeta(4);

  for (int q = 0; q < numruns; q++)
  {

    int currindexq = q * (nogen2 + 1) * lattsize * 6;
    aMeta.clear();
    // initialize lattice
    for (int i = 0; i < lattsize; i++)
    {
      for (int j = 0; j < lattsize2; j++)
      {
        latt[i][j][0] = 1;
        latt[i][j][1] = 0;
        latt[i][j][2] = 0;
        if (j > left && j <= right)
        {
          latt[i][j][0] = 2;
        }
        int m = latt[i][j][0];
        if (j == left || j == left + 1 || j == right || j == right + 1)
        {
          active[aMeta[0]] = {i, j, m};
          latt[i][j][1] = 1;
          latt[i][j][2] = aMeta[0];
          aMeta[0]++;
          aMeta[m]++;
        }
      }
    }

    t2 = 0;
    double t = 0.0;
    while (t <= pow(2, nogen2) + 1)
    {

      if (t > 0.0 && log2(t) >= t2)
      {
        pow_2 = true;
      }
      else
      {
        pow_2 = false;
      }

      if (aMeta[0] != aMeta[1] + aMeta[2] + aMeta[3])
        cout << "Something wrong with aMeta: (" << aMeta[0] << "," << aMeta[1] << "," << aMeta[2] << "," << aMeta[3] << ")" << endl;

      //choose a random lattice site to update from the list of active sites
      rannum1 = dis(generator);
      active_rand = chooseactive(active, aMeta, rannum1, sG, sW);
      i_rand = active[active_rand][0];
      j_rand = active[active_rand][1];
      int a, b, c, d, e, f, g;
      int lxindex = mod(i_rand - 1, lattsize);
      int rxindex = mod(i_rand + 1, lattsize);
      int uyindex = mod(j_rand - 1, lattsize);
      int dyindex = mod(j_rand + 1, lattsize);
      if (j_rand % 2 == 0)
      {
        //find neighbors of the chosen site if j-index is even
        a = latt[lxindex][j_rand][0];
        b = latt[rxindex][j_rand][0];
        c = latt[i_rand][uyindex][0];
        d = latt[rxindex][uyindex][0];
        e = latt[i_rand][dyindex][0];
        f = latt[rxindex][dyindex][0];
        g = latt[i_rand][j_rand][0];
      }
      else
      {
        //find neighbors of the chosen site if j-index is odd
        a = latt[lxindex][j_rand][0];
        b = latt[rxindex][j_rand][0];
        c = latt[i_rand][uyindex][0];
        d = latt[lxindex][uyindex][0];
        e = latt[i_rand][dyindex][0];
        f = latt[lxindex][dyindex][0];
        g = latt[i_rand][j_rand][0];
      }
      int winner = 0;
      if (a == b && b == c && c == d && d == e && e == f)
      {
        //if all the neighbors of the chosen site have the same identity, don't pull a random number; replace active site with identity of the neighbors
        if (a != g)
        {
          winner = a;
        }
        else
        {
          winner = g;
        }
      }
      else
      {
        //if the neighbors of the chosen cell have different identities, pull a random number and decide the winner to replace the active site
        rannum1 = dis(generator);
        winner = winstrain(a, b, c, d, e, f, gammas, rannum1);
      }
      if (winner == 2 && mu > 0.0000000001)
      {
        //if mutation rate is non-zero, check for a mutation event
        rannum2 = dis(generator);
        if (rannum2 < mu)
          winner = 3;
      }
      if (winner != g)
      {
        //if the lattice has changed, update the lattice and active list meta data
        latt[i_rand][j_rand][0] = winner;
        active[active_rand][2] = winner;
        --aMeta[g];
        ++aMeta[winner];
        //if the lattice has changed during the last update, update active list
        vector<vector<int>> tmp(7, vector<int>(3));
        if (j_rand % 2 == 0)
        {
          //find neighbors if j-index for current updated cell is even
          tmp = {{lxindex, j_rand, a},
                 {rxindex, j_rand, b},
                 {i_rand, uyindex, c},
                 {rxindex, uyindex, d},
                 {i_rand, dyindex, e},
                 {rxindex, dyindex, f},
                 {i_rand, j_rand, latt[i_rand][j_rand][0]}};
        }
        else
        {
          //find neighbors if j-index for current updated cell is odd
          tmp = {{lxindex, j_rand, a},
                 {rxindex, j_rand, b},
                 {i_rand, uyindex, c},
                 {lxindex, uyindex, d},
                 {i_rand, dyindex, e},
                 {lxindex, dyindex, f},
                 {i_rand, j_rand, latt[i_rand][j_rand][0]}};
        }
        for (auto &vec : tmp)
        {
          //for each of the neighbors of the updated cell (including updated cell), check if that cell is now active or inactive
          int new_i = vec[0], new_j = vec[1], new_m = vec[2];
          int in_active = latt[new_i][new_j][1];
          int index = latt[new_i][new_j][2];
          int ta, tb, tc, td, te, tf;
          int n_lxindex = mod(new_i - 1, lattsize);
          int n_rxindex = mod(new_i + 1, lattsize);
          int n_uyindex = mod(new_j - 1, lattsize);
          int n_dyindex = mod(new_j + 1, lattsize);
          if (new_j % 2 == 0)
          {
            //find nearest neighbors of site (new_i,new_j) if new_j is even
            ta = latt[n_lxindex][new_j][0];
            tb = latt[n_rxindex][new_j][0];
            tc = latt[new_i][n_uyindex][0];
            td = latt[n_rxindex][n_uyindex][0];
            te = latt[new_i][n_dyindex][0];
            tf = latt[n_rxindex][n_dyindex][0];
          }
          else
          {
            //find nearest neighbors of site (new_i,new_j) if new_j is odd
            ta = latt[n_lxindex][new_j][0];
            tb = latt[n_rxindex][new_j][0];
            tc = latt[new_i][n_uyindex][0];
            td = latt[n_lxindex][n_uyindex][0];
            te = latt[new_i][n_dyindex][0];
            tf = latt[n_lxindex][n_dyindex][0];
          }
          //check if current site is in the active list
          if (in_active == 1)
          {
            //if current site is in the active list, check if it *should* be in the active list
            if (new_m == ta && new_m == tb && new_m == tc && new_m == td && new_m == te && new_m == tf)
            {
              //if current site should not be in active list anymore (has same identity as all its neihbors), remove it
              active[index] = active[aMeta[0] - 1];
              latt[active[index][0]][active[index][1]][2] = index;
              latt[new_i][new_j][1] = 0;
              latt[new_i][new_j][2] = active.size() + 1;
              --aMeta[0];
              --aMeta[new_m];
            }
          }
          else if (in_active == 0)
          {
            if (new_m != ta || new_m != tb || new_m != tc || new_m != td || new_m != te || new_m != tf)
            {
              //if current lattice site is not in the active list and has at least one neighbor with a different identity, add it to active list
              active[aMeta[0]] = {new_i, new_j, new_m};
              latt[new_i][new_j][1] = 1;
              latt[new_i][new_j][2] = aMeta[0];
              ++aMeta[0];
              ++aMeta[new_m];
            }
          }
        }
      }

      if (pow_2)
      {

        int currindexqt2 = currindexq + t2 * lattsize * 6;
        cout << ": (q,t2,t) = (" << q << "," << t2 << "," << t << "); aSize = " << aMeta[0] << endl;

        for (int i = 0; i < lattsize; i++)
        {

          //Find first and last elements in current column of the lattice with value >1
          for (auto it = latt[i].begin(); it != latt[i].end(); ++it)
          {
            if ((*it)[0] > 1)
            {
              Lit = it;
              break;
            }
          }
          for (auto it = latt[i].rbegin(); it != latt[i].rend(); it++)
          {
            if ((*it)[0] > 1)
            {
              Rit = it;
              break;
            }
          }
          //The first and last elements in the current column with value >1 will act as our maximum range for this column
          //since the the identity for j<Lmax and j>Ramax must be 1
          Lmax = distance(latt[i].begin(), Lit);
          Rmax = distance(latt[i].begin(), (Rit + 1).base());

          std::vector<int> Lmatches;
          std::vector<int> Rmatches;
          //Find all green cells (M = 1) within the range Lmax < j < Rmax and add these cells to the vector Lmatches or Rmatches
          //depending on which side of the lattice the green cell is on
          for (auto j = Lit, toofar = (Rit + 1).base() + 1; j != toofar; ++j)
          {
            if ((*j)[0] == 1 && distance(latt[i].begin(), j) < lattsize2 / 2)
              Lmatches.push_back(distance(latt[i].begin(), j));
            if ((*j)[0] == 1 && distance(latt[i].begin(), j) > lattsize2 / 2)
              Rmatches.push_back(distance(latt[i].begin(), j));
          }

          double Lavg = 0, Ravg = 0;
          //average the j-indices of the green cells for each column to determine the average "position" of the boundary in this column
          if (Lmatches.size() > 0)
          {
            LInt = Lmatches[Lmatches.size() - 1] - Lmax;
            Lavg = (Lmax - 1.0) / ((double)Lmatches.size() + 1.0);
            for (auto k = Lmatches.begin(); k != Lmatches.end(); ++k)
            {
              Lavg += *k / ((double)Lmatches.size() + 1.0);
            }
            if (Lavg < 0.000001)
              Lavg = Lmax;
          }
          else
          {
            Lavg = Lmax;
            LInt = 0;
          }
          if (Rmatches.size() > 0)
          {
            RInt = Rmax - Rmatches[0];
            Ravg = (Rmax + 1.0) / ((double)Rmatches.size() + 1.0);
            for (auto k = Rmatches.begin(); k != Rmatches.end(); ++k)
            {
              Ravg += *k / ((double)Rmatches.size() + 1.0);
            }
            if (Ravg < 0.000001)
              Ravg = Rmax;
          }
          else
          {
            Ravg = Rmax;
            RInt = 0;
          }
          int currindex = currindexqt2 + i * 6;
          results[currindex + 0] = Lavg;
          results[currindex + 1] = Ravg;
          results[currindex + 2] = Lmax;
          results[currindex + 3] = Rmax;
          results[currindex + 4] = LInt;
          results[currindex + 5] = RInt;
        }
        pow_2 = false;
        t2++;
      }
      //Advance the time according to Gillespie's algorithm
      dt = -log(1.0 - dis(generator)) / (double)aMeta[0];
      t = t + dt;
    }
  }

  if (image)
  {
    for (auto it = latt.begin(); it != latt.end(); ++it)
    {
      for (auto jt = (*it).begin(); jt != (*it).end(); ++jt)
      {
        double x = (double)distance(latt.begin(), it) - 0.5 * mod(distance((*it).begin(), jt), 2);
        double y = (double)distance((*it).begin(), jt) * sqrt(3) / 2;
        if ((*jt)[0] > 1)
          outstats3 << t2 << "\t" << x << "\t" << y << "\t" << (*jt)[0] << "\t" << (*jt)[1] << endl;
      }
    }
    outstats3.close();
  }

  //delete latt and active vectors to clear memory for new vectors for analysis
  latt.clear();
  vector<vector<vector<int>>>().swap(latt);
  active.clear();
  vector<vector<int>>().swap(active);

  vector<vector<vector<double>>> avgs(numruns, vector<vector<double>>(nogen2 + 1, vector<double>(6)));
  vector<vector<vector<double>>> avgs2(numruns, vector<vector<double>>(nogen2 + 1, vector<double>(6)));

  for (int i = 0; i < numruns; i++)
  {
    int currindexq = i * (nogen2 + 1) * lattsize * 6;
    for (int j = 0; j <= nogen2; j++)
    {
      int currindexqt2 = currindexq + j * lattsize * 6;
      for (int k = 0; k < lattsize; k++)
      {
        int currindex = currindexqt2 + k * 6;
        for (int l = 0; l < 6; l++)
          avgs[i][j][l] += results[currindex + l] / (double)lattsize;
      }
    }
  }
  for (int i = 0; i < numruns; i++)
  {
    int currindexq = i * (nogen2 + 1) * lattsize * 6;
    for (int j = 0; j <= nogen2; j++)
    {
      int currindexqt2 = currindexq + j * lattsize * 6;
      for (int k = 0; k < lattsize; k++)
      {
        int currindex = currindexqt2 + k * 6;
        for (int l = 0; l < 6; l++)
          avgs2[i][j][l] += pow(avgs[i][j][l] - results[currindex + l], 2) / (double)lattsize;
      }
    }
  }

  //clear results vector to clear memory for new vectors
  results.clear();
  vector<double>().swap(results);

  vector<double> w_avgs((nogen2 + 1) * 3);
  vector<double> sw_avgs((nogen2 + 1) * 3);

  for (int i = 0; i < numruns; i++)
  {
    for (int j = 0; j <= nogen2; j++)
    {
      for (int k = 0; k < 3; k++)
        w_avgs[j * 3 + k] += (avgs2[i][j][2 * k] + avgs2[i][j][2 * k + 1]) / 2.0 / (double)numruns;
    }
  }

  //only print results vector for run q = 0, t2 = nogen2
  int currindext2 = nogen2 * lattsize * 6;
  for (int k = 0; k < lattsize; k++)
  {
    int currindex = currindext2 + k * 6;
    outstats << "0\t" << nogen2 << "\t" << k << "\t" << results[currindex + 0]
             << "\t" << results[currindex + 1]
             << "\t" << results[currindex + 2]
             << "\t" << results[currindex + 3]
             << "\t" << results[currindex + 4]
             << "\t" << results[currindex + 5] << endl;
  }
  for (int i = 0; i <= nogen2; i++)
    outstats2 << pow(2, i) << "\t" << sqrt(sw_avgs[i * 3 + 0])
              << "\t" << sqrt(sw_avgs[i * 3 + 1])
              << "\t" << sqrt(sw_avgs[i * 3 + 2]) << endl;

  outstats.close();

  outstats2.close();

  return 0;
}
