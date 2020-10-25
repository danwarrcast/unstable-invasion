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
int winstrain(int a, int b, int c, int d, int e, int f, vector<double> &gammas, double rannum)
{

  double G = 0.0;
  int index = (a - 1) * pow(3, 6) + (b - 1) * pow(3, 5) + (c - 1) * pow(3, 4) + (d - 1) * pow(3, 3) + (e - 1) * pow(3, 2) + (f - 1) * 3;
  int winner;

  for (int s = 0; s < 3; s++)
  {
    G += gammas[index + s];
    if (G >= rannum)
    {
      winner = s + 1;
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

vector<int> chooseactive(vector<int> &aM, double rN1, double rN2, double sG, double sW)
{
  int winner;
  double sum = 0;
  double denom = (1.0 - sG) * (double)aM[1] + (double)aM[2] + (1.0 - sW) * (double)aM[3];
  vector<double> sums = {0, (1.0 - sG) * (double)aM[1] / denom, (double)aM[2] / denom, (1.0 - sW) * (double)aM[3] / denom};
  for (int i = 1; i < 4; i++)
  {
    sum += sums[i];
    if (sum >= rN1)
    {
      winner = i;
      break;
    }
  }
  int winner_index = round(rN2 * (aM[winner] - 1));
  return {winner_index, winner};
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
  vector<double> gammas(3 * 3 * 3 * 3 * 3 * 3 * 3);
  for (int i = 0; i < 3; i++)
  {
    int curri = i * 3 * 3 * 3 * 3 * 3 * 3;
    for (int j = 0; j < 3; j++)
    {
      int currj = curri + j * 3 * 3 * 3 * 3 * 3;
      for (int k = 0; k < 3; k++)
      {
        int currk = currj + k * 3 * 3 * 3 * 3;
        for (int l = 0; l < 3; l++)
        {
          int currl = currk + l * 3 * 3 * 3;
          for (int m = 0; m < 3; m++)
          {
            int currm = currl + m * 3 * 3;
            for (int n = 0; n < 3; n++)
            {
              vector<int> tmp = {0, 0, 0};
              tmp[i]++;
              tmp[j]++;
              tmp[k]++;
              tmp[l]++;
              tmp[m]++;
              tmp[n]++;
              double Denom = (1.0 - sG) * (double)tmp[0] + (double)tmp[1] + (1.0 - sW) * (double)tmp[2];
              gammas[currm + n * 3 + 0] = (1.0 - sG) * (double)tmp[0] / Denom;
              gammas[currm + n * 3 + 1] = (double)tmp[1] / Denom;
              gammas[currm + n * 3 + 2] = (1.0 - sW) * (double)tmp[2] / Denom;
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
  int Lmax, Rmax, LInt, RInt, i_rand, j_rand;
  vector<int> active_rand;
  double dt;

  char *_emergencyMemory = new char[16384];
  try
  {
    vector<double> results(numruns * (nogen2 + 1) * lattsize * 6);
    vector<int> latt(lattsize * lattsize2 * 3);
    vector<int> active(3 * lattsize * lattsize2 * 2);
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
  vector<int> latt(lattsize * lattsize2 * 3);
  vector<int> active(3 * lattsize * lattsize2 * 2);
  vector<int> aMeta(4);

  for (int q = 0; q < numruns; q++)
  {

    int currindexq = q * (nogen2 + 1) * lattsize * 6;
    fill(aMeta.begin(), aMeta.end(), 0);
    // initialize lattice
    int atmp = lattsize * lattsize2 * 2;
    for (int i = 0; i < lattsize; i++)
    {
      int tmpindex = i * lattsize2 * 3;
      for (int j = 0; j < lattsize2; j++)
      {
        latt[tmpindex + j * 3 + 0] = 1;
        latt[tmpindex + j * 3 + 1] = 0;
        latt[tmpindex + j * 3 + 2] = 0;
        if (j > left && j <= right)
        {
          latt[tmpindex + j * 3 + 0] = 2;
        }
        int m = latt[tmpindex + j * 3 + 0];
        if (j == left || j == left + 1 || j == right || j == right + 1)
        {
          active[(m - 1) * atmp + aMeta[m] * 2] = i;
          active[(m - 1) * atmp + aMeta[m] * 2 + 1] = j;
          latt[tmpindex + j * 3 + 1] = 1;
          latt[tmpindex + j * 3 + 2] = aMeta[m];
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
      rannum2 = dis(generator);
      //active_rand = round(rannum1*(aMeta[0]-1));
      active_rand = chooseactive(aMeta, rannum1, rannum2, sG, sW);
      int chosenIndex = active_rand[0];
      int g = active_rand[1];
      i_rand = active[(g-1) * atmp + chosenIndex * 2];
      j_rand = active[(g-1) * atmp + chosenIndex * 2 + 1];
      int a, b, c, d, e, f;
      int lxindex = mod(i_rand - 1, lattsize);
      int rxindex = mod(i_rand + 1, lattsize);
      int uyindex = mod(j_rand + 1, lattsize2);
      int dyindex = mod(j_rand - 1, lattsize2);
      if (j_rand % 2 == 0)
      {
        //find neighbors of the chosen site if j-index is even
        a = latt[lxindex * lattsize2 * 3 + j_rand * 3];
        b = latt[rxindex * lattsize2 * 3 + j_rand * 3];
        c = latt[i_rand * lattsize2 * 3 + uyindex * 3];
        d = latt[rxindex * lattsize2 * 3 + uyindex * 3];
        e = latt[i_rand * lattsize2 * 3 + dyindex * 3];
        f = latt[rxindex * lattsize2 * 3 + dyindex * 3];
      }
      else
      {
        //find neighbors of the chosen site if j-index is odd
        a = latt[lxindex * lattsize2 * 3 + j_rand * 3];
        b = latt[rxindex * lattsize2 * 3 + j_rand * 3];
        c = latt[i_rand * lattsize2 * 3 + uyindex * 3];
        d = latt[lxindex * lattsize2 * 3 + uyindex * 3];
        e = latt[i_rand * lattsize2 * 3 + dyindex * 3];
        f = latt[lxindex * lattsize2 * 3 + dyindex * 3];
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
        latt[i_rand * lattsize2 * 3 + j_rand * 3] = winner;
        active[(winner-1) * atmp + aMeta[winner] * 2] = i_rand;
        active[(winner-1) * atmp + aMeta[winner] * 2 + 1] = j_rand;
        latt[i_rand * lattsize2 * 3 + j_rand * 3 + 2] = aMeta[winner];
        active[(g-1) * atmp + chosenIndex * 2] = active[(g-1) * atmp + aMeta[g] * 2];
        active[(g-1) * atmp + chosenIndex * 2 + 1] = active[(g-1) * atmp + aMeta[g] * 2 + 1];
        latt[active[(g-1) * atmp + chosenIndex * 2] * lattsize2 * 3 + active[(g-1) * atmp + chosenIndex * 2 + 1] * 3 + 2] = chosenIndex;
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
                 {i_rand, j_rand, winner}};
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
                 {i_rand, j_rand, winner}};
        }
        for (auto &vec : tmp)
        {
          //for each of the neighbors of the updated cell (including updated cell), check if that cell is now active or inactive
          int new_i = vec[0], new_j = vec[1], new_m = vec[2];
          int in_active = latt[new_i * lattsize2 * 3 + new_j * 3 + 1];
          int index = latt[new_i * lattsize2 * 3 + new_j * 3 + 2];
          int ta, tb, tc, td, te, tf;
          int n_lxindex = mod(new_i - 1, lattsize);
          int n_rxindex = mod(new_i + 1, lattsize);
          int n_uyindex = mod(new_j + 1, lattsize2);
          int n_dyindex = mod(new_j - 1, lattsize2);
          if (new_j % 2 == 0)
          {
            //find nearest neighbors of site (new_i,new_j) if new_j is even
            ta = latt[n_lxindex * lattsize2 * 3 + new_j * 3];
            tb = latt[n_rxindex * lattsize2 * 3 + new_j * 3];
            tc = latt[new_i * lattsize2 * 3 + n_uyindex * 3];
            td = latt[n_rxindex * lattsize2 * 3 + n_uyindex * 3];
            te = latt[new_i * lattsize2 * 3 + n_dyindex * 3];
            tf = latt[n_rxindex * lattsize2 * 3 + n_dyindex * 3];
          }
          else
          {
            //find nearest neighbors of site (new_i,new_j) if new_j is odd
            ta = latt[n_lxindex * lattsize2 * 3 + new_j * 3];
            tb = latt[n_rxindex * lattsize2 * 3 + new_j * 3];
            tc = latt[new_i * lattsize2 * 3 + n_uyindex * 3];
            td = latt[n_lxindex * lattsize2 * 3 + n_uyindex * 3];
            te = latt[new_i * lattsize2 * 3 + n_dyindex * 3];
            tf = latt[n_lxindex * lattsize2 * 3 + n_dyindex * 3];
          }
          //check if current site is in the active list
          if (in_active == 1)
          {
            //if current site is in the active list, check if it *should* be in the active list
            if (new_m == ta && new_m == tb && new_m == tc && new_m == td && new_m == te && new_m == tf)
            {
              //if current site should not be in active list anymore (has same identity as all its neihbors), remove it
              active[(new_m-1) * atmp + index * 2] = active[(new_m-1) * atmp + (aMeta[new_m]-1) * 2];
              active[(new_m-1) * atmp + index * 2 + 1] = active[(new_m-1) * atmp + (aMeta[new_m]-1) * 2 + 1];
              latt[active[(new_m-1) * atmp + index * 2] * lattsize2 * 3 + active[(new_m-1) * atmp + index * 2 + 1] * 3 + 2] = index;
              latt[new_i * lattsize2 * 3 + new_j * 3 + 1] = 0;
              latt[new_i * lattsize2 * 3 + new_j * 3 + 2] = atmp;
              --aMeta[0];
              --aMeta[new_m];
            }
          }
          else if (in_active == 0)
          {
            if (new_m != ta || new_m != tb || new_m != tc || new_m != td || new_m != te || new_m != tf)
            {
              //if current lattice site is not in the active list and has at least one neighbor with a different identity, add it to active list
              active[(new_m-1) * atmp + aMeta[new_m] * 2] = new_i;
              active[(new_m-1) * atmp + aMeta[new_m] * 2 + 1] = new_j;
              latt[new_i * lattsize2 * 3 + new_j * 3 + 1] = 1;
              latt[new_i * lattsize2 * 3 + new_j * 3 + 2] = aMeta[new_m];
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

        for (int i = 0; i < lattsize; ++i)
        {
          int tmpindex = i * lattsize2 * 3;
          //Find first and last elements in current column of the lattice with value >1
          for (int j = 0; j < lattsize2; ++j)
          {
            if (latt[tmpindex + j * 3] > 1)
            {
              Lmax = j;
              break;
            }
          }
          for (int j = lattsize2 - 1; j >= 0; --j)
          {
            if (latt[tmpindex + j * 3] > 1)
            {
              Rmax = j;
              break;
            }
          }

          std::vector<int> Lmatches;
          std::vector<int> Rmatches;
          //Find all green cells (M = 1) within the range Lmax < j < Rmax and add these cells to the vector Lmatches or Rmatches
          //depending on which side of the lattice the green cell is on
          for (int j = Lmax; j < Rmax; j++)
          {
            if (latt[tmpindex + j * 3] == 1 && j < lattsize2 / 2)
              Lmatches.push_back(j);
            if (latt[tmpindex + j * 3] == 1 && j > lattsize2 / 2)
              Rmatches.push_back(j);
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
    for (int i = 0; i < lattsize; ++i)
    {
      int tmpindex = i * lattsize2 * 3;
      for (int j = 0; j < lattsize2; ++j)
      {
        double x = (double)i - 0.5 * (j % 2);
        double y = j * sqrt(3) / 2;
        if (latt[tmpindex + j * 3] > 1)
        {
          outstats3 << t2 << "\t" << x << "\t" << y << "\t" << latt[tmpindex + j * 3] << endl;
        }
      }
    }
    outstats3.close();
  }

  //delete latt and active vectors to clear memory for new vectors for analysis
  latt.clear();
  vector<int>().swap(latt);
  active.clear();
  vector<int>().swap(active);

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
    outstats.close();

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

    for (int i = 0; i <= nogen2; i++)
      outstats2 << pow(2, i) << "\t" << sqrt(sw_avgs[i * 3 + 0])
                << "\t" << sqrt(sw_avgs[i * 3 + 1])
                << "\t" << sqrt(sw_avgs[i * 3 + 2]) << endl;


  return 0;
}
