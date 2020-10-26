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
  double denom = (1.0 + sG) * (double)aM[1] + (double)aM[2] + (1.0 + sW) * (double)aM[3];
  vector<double> sums = {0, (1.0 + sG) * (double)aM[1] / denom, (double)aM[2] / denom, (1.0 + sW) * (double)aM[3] / denom};
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
  outstats2.open(statsfilename.c_str());

  cout << statsfilename.c_str() << endl;

  statsfilename = "diffusion_out/man/ranseq_boundary";

  tempstring.str("");
  tempstring.clear();
  tempstring << "_run" << filecounter;
  statsfilename += tempstring.str();

  outstats.open(statsfilename.c_str());

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

  outstats2 << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " Has_Image: " << image << endl;
  outstats2 << " # sW = " << sW << " b = " << sW - sG << " mu = " << mu << endl;

  if (image)
  {
    outstats3 << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " Has_Image: " << image << endl;
    outstats3 << " # sW = " << sW << " b = " << sW - sG << " mu = " << mu << endl;
  }

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
  int nogen2 = ceil(log2(nogen));
  int lattsize2 = lattsize * 4;

  int left = (lattsize2) / 2 - lattsize / 4, right = (lattsize2) / 2 + lattsize / 4;
  int Lmax, Rmax, LInt, RInt;

  char *_emergencyMemory = new char[16384];
  try
  {
    vector<double> results(numruns * (nogen2 + 1) * lattsize * lattsize * 6);
    vector<int> latt(lattsize * lattsize * lattsize2 * 3);
    vector<int> active(3 * lattsize * lattsize * lattsize2 * 3);
    vector<int> aMeta(4);
  }
  catch (bad_alloc &ex)
  {
    delete[] _emergencyMemory;
    cerr << "Not enough memory!" << endl;
    cin.get();
    exit(1);
  }
  vector<double> results(numruns * (nogen2 + 1) * lattsize * lattsize * 6);
  vector<int> latt(lattsize * lattsize * lattsize2 * 3);
  vector<int> active(3 * lattsize * lattsize * lattsize2 * 3);
  vector<int> aMeta(4);

  for (int q = 0; q < numruns; q++)
  {

    int currindexq = q * (nogen2 + 1) * lattsize * lattsize * 6;
    fill(aMeta.begin(), aMeta.end(), 0);
    int atmp = lattsize * lattsize * lattsize2 * 3;

    // initialize lattice
    for (int i = 0; i < lattsize; i++)
    {
      int tmpindex = i * lattsize * lattsize2 * 3;
      for (int j = 0; j < lattsize; j++)
      {
        int tmpindex2 = tmpindex + j * lattsize2 * 3;
        for (int k = 0; k < lattsize2; k++)
        {
          int tmpindex3 = tmpindex2 + k * 3;
          latt[tmpindex3] = 1;
          latt[tmpindex3 + 1] = 0;
          latt[tmpindex3 + 2] = 0;
          if (k > left && k <= right)
          {
            latt[tmpindex3] = 2;
          }
          int m = latt[tmpindex2];
          if (k == left || k == left + 1 || k == right || k == right + 1)
          {
            int tmpindexa = (m - 1) * atmp + aMeta[m] * 3;
            active[tmpindexa] = i;
            active[tmpindexa + 1] = j;
            active[tmpindexa + 2] = k;
            latt[tmpindex3 + 1] = 1;
            latt[tmpindex3 + 2] = aMeta[m];
            aMeta[0]++;
            aMeta[m]++;
          }
        }
      }
    }

    int t2 = 0;
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
        cout << "t = " << t << ": Something wrong with aMeta: (" << aMeta[0] << "," << aMeta[1] << "," << aMeta[2] << "," << aMeta[3] << ")" << endl;

      //choose a random lattice site to update from the list of active sites
      rannum1 = dis(generator);
      rannum2 = dis(generator);
      //active_rand = round(rannum1*(aMeta[0]-1));
      vector<int> active_rand(2);
      active_rand = chooseactive(aMeta, rannum1, rannum2, sG, sW);
      int chosenIndex = active_rand[0];
      int g = active_rand[1];
      int tmpindexa = (g - 1) * atmp + chosenIndex * 3;
      int i_rand = active[tmpindexa];
      int j_rand = active[tmpindexa + 1];
      int k_rand = active[tmpindexa + 2];
      int a, b, c, d, e, f;
      int lxindex = mod(i_rand - 1, lattsize);
      int rxindex = mod(i_rand + 1, lattsize);
      int byindex = mod(j_rand - 1, lattsize);
      int fyindex = mod(j_rand + 1, lattsize);
      int uzindex = mod(k_rand + 1, lattsize2);
      int dzindex = mod(k_rand - 1, lattsize2);
      //find neighbors of the chosen site if j-index is even
      int lsxls2x3 = lsxls2x3;
      int ls2x3 = lattsize2 * 3;
      a = latt[lxindex * lsxls2x3 + j_rand * ls2x3 + k_rand * 3];
      b = latt[rxindex * lsxls2x3 + j_rand * ls2x3 + k_rand * 3];
      c = latt[i_rand * lsxls2x3 + byindex * ls2x3 + k_rand * 3];
      d = latt[i_rand * lsxls2x3 + fyindex * ls2x3 + k_rand * 3];
      e = latt[i_rand * lsxls2x3 + j_rand * ls2x3 + uzindex * 3];
      f = latt[i_rand * lsxls2x3 + j_rand * ls2x3 + dzindex * 3];
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

        //update lattice with new cell identity at the chosen site
        int tmpindex = i_rand * lsxls2x3 + j_rand * ls2x3 + k_rand * 3;
        latt[tmpindex] = winner;

        //add this site to the active list
        tmpindexa = (winner - 1) * atmp + aMeta[winner] * 3;
        active[tmpindexa] = i_rand;
        active[tmpindexa + 1] = j_rand;
        active[tmpindexa + 2] = k_rand;

        //tell lattice where to find its new location in the active lisr
        latt[tmpindex + 2] = aMeta[winner];

        //remove the updated site from its old active list
        tmpindexa = (g - 1) * atmp + chosenIndex * 3;
        int tmpindexa2 = (g - 1) * atmp + (aMeta[g] - 1) * 3;
        active[tmpindexa] = active[tmpindexa2];
        active[tmpindexa + 1] = active[tmpindexa2 + 1];
        active[tmpindexa + 2] = active[tmpindexa2 + 2];

        //tell the lattice about the changes we just made
        latt[active[tmpindexa] * lsxls2x3 + active[tmpindexa + 1] * ls2x3 + active[tmpindexa + 2] * 3 + 2] = chosenIndex;

        //update counts
        --aMeta[g];
        ++aMeta[winner];

        //if the lattice has changed during the last update, update active list
        vector<vector<int>> tmp(7, vector<int>(3));
        tmp = {{lxindex, j_rand, k_rand, a},
               {rxindex, j_rand, k_rand, b},
               {i_rand, byindex, k_rand, c},
               {i_rand, fyindex, k_rand, d},
               {i_rand, j_rand, uzindex, e},
               {i_rand, j_rand, dzindex, f},
               {i_rand, j_rand, k_rand, winner}};
        for (auto &vec : tmp)
        {
          //for each of the neighbors of the updated cell (including updated cell), check if that cell is now active or inactive
          int new_i = vec[0], new_j = vec[1], new_k = vec[2], new_m = vec[3];
          tmpindex = new_i * lsxls2x3 + new_j * ls2x3 + new_k * 3;
          int in_active = latt[tmpindex + 1];
          int index = latt[tmpindex + 2];
          int ta, tb, tc, td, te, tf;
          int n_lxindex = mod(new_i - 1, lattsize);
          int n_rxindex = mod(new_i + 1, lattsize);
          int n_byindex = mod(new_j - 1, lattsize);
          int n_fyindex = mod(new_j + 1, lattsize);
          int n_uzindex = mod(new_k + 1, lattsize2);
          int n_dzindex = mod(new_k - 1, lattsize2);
          ta = latt[n_lxindex * lsxls2x3 + new_j * ls2x3 + new_k * 3];
          tb = latt[n_rxindex * lsxls2x3 + new_j * ls2x3 + new_k * 3];
          tc = latt[new_i * lsxls2x3 + n_byindex * ls2x3 + new_k * 3];
          td = latt[new_i * lsxls2x3 + n_fyindex * ls2x3 + new_k * 3];
          te = latt[new_i * lsxls2x3 + new_j * ls2x3 + uzindex * 3];
          tf = latt[new_i * lsxls2x3 + new_j * ls2x3 + dzindex * 3];
          //check if current site is in the active list
          if (in_active == 1)
          {
            //if current site is in the active list, check if it *should* be in the active list
            if (new_m == ta && new_m == tb && new_m == tc && new_m == td && new_m == te && new_m == tf)
            {
              //if current site should not be in active list anymore (has same identity as all its neihbors), remove it
              tmpindexa = (new_m - 1) * atmp + index * 3;
              tmpindexa2 = (new_m - 1) * atmp + (aMeta[new_m] - 1) * 3;
              active[tmpindexa] = active[tmpindexa2];
              active[tmpindexa + 1] = active[tmpindexa2 + 1];
              active[tmpindexa + 2] = active[tmpindexa2 + 2];
              latt[active[tmpindexa] * lsxls2x3 + active[tmpindexa + 1] * ls2x3 + active[tmpindexa + 2] * 3 + 2] = index;
              tmpindex = new_i * lsxls2x3 + new_j * ls2x3 + new_k * 3;
              latt[tmpindex + 1] = 0;
              latt[tmpindex + 2] = atmp;
              --aMeta[0];
              --aMeta[new_m];
            }
          }
          else if (in_active == 0)
          {
            if (new_m != ta || new_m != tb || new_m != tc || new_m != td || new_m != te || new_m != tf)
            {
              //if current lattice site is not in the active list and has at least one neighbor with a different identity, add it to active list
              tmpindexa = (new_m - 1) * atmp + aMeta[new_m] * 3;
              active[tmpindexa] = new_i;
              active[tmpindexa + 1] = new_j;
              active[tmpindexa + 2] = new_k;
              tmpindex = new_i * lsxls2x3 + new_j * ls2x3 + new_k * 3;
              latt[tmpindex + 1] = 1;
              latt[tmpindex + 2] = aMeta[new_m];
              ++aMeta[0];
              ++aMeta[new_m];
            }
          }
        }
      }

      if (pow_2)
      {

        int currindexqt2 = currindexq + t2 * lattsize * lattsize * 6;
        cout << ": (q,t2,t) = (" << q << "," << t2 << "," << t << "); aSize = " << aMeta[0] << endl;

        for (int i = 0; i < lattsize; ++i)
        {
          int tmpindex = i * lsxls2x3;
          int currindexi = currindexqt2 + i * lattsize * 6;
          //Find first and last elements in current column of the lattice with value >1
          for (int j = 0; j < lattsize; ++j)
          {
            int tmpindex2 = tmpindex + j * ls2x3;
            int currindex = currindexi + j * 6;
            for (int k = 0; k < lattsize2; ++k)
            {
              if (latt[tmpindex2 + k * 3] > 1)
              {
                Lmax = j;
                break;
              }
            }
            for (int k = lattsize2 - 1; k >= 0; --k)
            {
              if (latt[tmpindex2 + k * 3] > 1)
              {
                Rmax = j;
                break;
              }
            }

            std::vector<int> Lmatches;
            std::vector<int> Rmatches;
            //Find all green cells (M = 1) within the range Lmax < j < Rmax and add these cells to the vector Lmatches or Rmatches
            //depending on which side of the lattice the green cell is on
            for (int l = Lmax; l < Rmax; ++l)
            {
              if (latt[tmpindex2 + l * 3] == 1 && l < lattsize2 / 2)
                Lmatches.push_back(l);
              if (latt[tmpindex2 + l * 3] == 1 && l > lattsize2 / 2)
                Rmatches.push_back(l);
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
            results[currindex + 0] = Lavg;
            results[currindex + 1] = Ravg;
            results[currindex + 2] = Lmax;
            results[currindex + 3] = Rmax;
            results[currindex + 4] = LInt;
            results[currindex + 5] = RInt;
          }
        }
        pow_2 = false;
        t2++;
      }
      //Advance the time according to Gillespie's algorithm
      double dt = -log(1.0 - dis(generator)) / (double)aMeta[0];
      t = t + dt;
    }
  }

  if (image)
  {
    for (int i = 0; i < lattsize; ++i)
    {
      int tmpindex = i * lattsize * lattsize2 * 3;
      for (int j = 0; j < lattsize; ++j)
      {
        int tmpindex2 = tmpindex + j * lattsize2 * 3;
        for (int k = 0; k < lattsize2; ++k)
        {
          int currindex = tmpindex2 + k * 3;
          if (latt[currindex] > 1)
          {
            outstats3 << nogen2 << "," << i << "," << j << "," << k << "," << latt[currindex] << endl;
          }
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

  vector<double> avgs(numruns * (nogen2 + 1) * 6);
  vector<double> avgs2(numruns * (nogen2 + 1) * 6);

  int ng2x6 = (nogen2 + 1) * 6;
  for (int i = 0; i < numruns; i++)
  {
    int currindexq = i * (nogen2 + 1) * lattsize * lattsize * 6;
    for (int j = 0; j <= nogen2; j++)
    {
      int currindexqt2 = currindexq + j * lattsize * lattsize * 6;
      for (int k = 0; k < lattsize; k++)
      {
        int currindexk = currindexqt2 + k * lattsize * 6;
        for (int l = 0; l < lattsize; l++)
        {
          int currindex = currindexk + l * 6;
          for (int s = 0; s < 6; ++s)
          {
            avgs[i * ng2x6 + j * 6 + s] += results[currindex + s] / (double)lattsize / (double)lattsize;
          }
        }
      }
    }
  }
  for (int i = 0; i < numruns; i++)
  {
    int currindexq = i * (nogen2 + 1) * lattsize * lattsize * 6;
    for (int j = 0; j <= nogen2; j++)
    {
      int currindexqt2 = currindexq + j * lattsize * lattsize * 6;
      for (int k = 0; k < lattsize; k++)
      {
        int currindexk = currindexqt2 + k * lattsize * 6;
        for (int l = 0; l < lattsize; l++)
        {
          int currindex = currindexk + l * 6;
          for (int s = 0; s < 6; ++s)
          {
            avgs2[i * ng2x6 + j * 6 + s] += pow(avgs[i * ng2x6 + j * 6 + s] - results[currindex + s], 2) / (double)lattsize / (double)lattsize;
          }
        }
      }
    }
  }

  avgs.clear();
  vector<double>().swap(avgs);

  //only print results vector for run q = 0, t2 = nogen2
  int currindext2 = nogen2 * lattsize * lattsize * 6;
  for (int k = 0; k < lattsize; k++)
  {
    int currindexk = currindext2 + k * lattsize * 6;
    for (int l = 0; l < lattsize; l++)
    {
      int currindex = currindexk + l * 6;
      outstats << "0," << nogen2 << "," << k << "," << results[currindex]
               << "," << results[currindex + 1]
               << "," << results[currindex + 2]
               << "," << results[currindex + 3]
               << "," << results[currindex + 4]
               << "," << results[currindex + 5] << endl;
    }
  }
  outstats.close();

  //clear results vector to clear memory for new vectors
  results.clear();
  vector<double>().swap(results);

  vector<double> w_avgs((nogen2 + 1) * 3);

  for (int i = 0; i < numruns; i++)
  {
    for (int j = 0; j <= nogen2; j++)
    {
      for (int k = 0; k < 3; k++)
        w_avgs[j * 3 + k] += (avgs2[i * ng2x6 * j * 6 + (2 * k)] + avgs2[i * ng2x6 * j * 6 + (2 * k + 1)]) / 2.0 / (double)numruns;
    }
  }
  for (int i = 0; i <= nogen2; i++)
  {
    outstats2 << pow(2, i) << "," << sqrt(w_avgs[i * 3 + 0])
              << "," << sqrt(w_avgs[i * 3 + 1])
              << "," << sqrt(w_avgs[i * 3 + 2]) << endl;
  }

  return 0;
}
