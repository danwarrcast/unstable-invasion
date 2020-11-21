#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <curand_kernel.h>
#include <algorithm>

bool to_bool(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

__global__
void setup_kernel(curandState *state, int N, unsigned long SEED) 
{        
    int i, idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    // Init Random vector
    for (i = idx; i < N; i += stride)
    {
        curand_init(SEED + i, 0, 0, &state[i]);
    }

}

__device__
int winstrain(int a, int b, int c, double s1, double s3, double rannum)
{
    double g_a = 1.0 - s3*(double)(a-1)*(double)(a-2)/2.0 - s1*(double)(2-a)*(double)(3-a)/2.0;
    double g_b = 1.0 - s3*(double)(b-1)*(double)(b-2)/2.0 - s1*(double)(2-b)*(double)(3-b)/2.0;
    double g_c = 1.0 - s3*(double)(c-1)*(double)(c-2)/2.0 - s1*(double)(2-c)*(double)(3-c)/2.0;
  
    double vec[6] = {g_a/(g_a+g_b+g_c), double(a), g_b/(g_a+g_b+g_c), double(b), g_c/(g_a+g_b+g_c), double(c)};
  
    double G = 0.0;
    double result;
  
    for (int s = 0; s < 3; s++) {
      G += vec[2*s];
      if (G > rannum) {
        result = (int)vec[2*s + 1];
        break;
      }
    }
  
    return result;
}

__global__
void update_odd(int N, int *l_d, int *l_u, int L, int L2, double s1, double s3, double mu, curandState *state)
{
    curandState localState;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += stride)
    {
        int x = index % L;
        int y = index / L;

        int yp = y<L2-1?y+1:0;
        //int ym = y>0?y-1:L2-1; not used in odd updates
        int xp = x<L-1?x+1:0;
        int xm = x>0?x-1:L-1;

        if (y % 2 == 0)
        {
           int a = l_d[index];
           int b = l_d[yp * L + x];
           int c = l_d[yp * L + xp];
           if (a == b && a == c)
           {
                l_u[index] = a;
           } else {
               localState = state[i];
               l_u[index] = winstrain(a, b, c, s1, s3, curand_uniform(&localState));
               state[i] = localState;
           }
        } else {
           int a = l_d[index];
           int b = l_d[yp * L + x];
           int c = l_d[yp * L + xm];
           if (a == b && a == c)
           {
                l_u[index] = a;
           } else {
               localState = state[i];
               l_u[index] = winstrain(a, b, c, s1, s3, curand_uniform(&localState));
               state[i] = localState;
           }
        }
        //if mu is non-zero, check for mutation event
        if (l_u[index] == 2 && mu > 0.0000000001)
        {
            localState = state[i];
            if (curand_uniform(&localState) < mu) l_u[index] = 3;
            state[i] = localState;
        }
    }
}

__global__
void update_even(int N, int *l_d, int *l_u, int L, int L2, double s1, double s3, double mu, curandState *state)
{
    curandState localState;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += stride)
    {
        int x = index % L;
        int y = index / L;

        //int yp = y<L2-1?y+1:0; not used in even updates
        int ym = y>0?y-1:L2-1;
        int xp = x<L-1?x+1:0;
        int xm = x>0?x-1:L-1;

        if (y % 2 == 0)
        {
           int a = l_d[index];
           int b = l_d[ym * L + x];
           int c = l_d[ym * L + xp];
           if (a == b && a == c)
           {
                l_u[index] = a;
           } else {
               localState = state[i];
               l_u[index] = winstrain(a, b, c, s1, s3, curand_uniform(&localState));
               state[i] = localState;
           }
        } else {
           int a = l_d[index];
           int b = l_d[ym * L + x];
           int c = l_d[ym * L + xm];
           if (a == b && a == c)
           {
                l_u[index] = a;
           } else {
               localState = state[i];
               l_u[index] = winstrain(a, b, c, s1, s3, curand_uniform(&localState));
               state[i] = localState;
           }
        }
        //if mu is non-zero, check for mutation event
        if (l_u[index] == 2 && mu > 0.0000000001)
        {
            localState = state[i];
            if (curand_uniform(&localState) < mu) l_u[index] = 3;
            state[i] = localState;
        }
    }
}

int main(int argc, char* argv[])
{
    cudaError_t cudaStatus;

    int lattsize, lattsize2, nogen, numruns;
    double mu, s1, s3;
    bool image;

    lattsize = atoi(argv[1]);
    lattsize2 = lattsize * 4;
    s3 = atof(argv[2]);
    s1 = s3 - atof(argv[3]);
    mu = atof(argv[4]);
    nogen = atoi(argv[5]);
    numruns = atoi(argv[6]);
    image = false;

    int N = lattsize * lattsize2;
    int *lattdown;
    int *lattup;
    curandState* devStates;

    int left = lattsize2 / 2 - lattsize / 4;
    int right = lattsize2 / 2 + lattsize / 4;

    int nogen2 = floor(log2(nogen));
    nogen = pow(2, nogen2);

    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    if (argc == 8 && to_bool(argv[7]))
    {
      image = true;
    }
    else if (argc == 8 && !to_bool(argv[7]))
    {
      image = false;
    }

    std::string statsfilename;  
    std::ofstream outstats;
    std::ofstream outstats2;
    std::ofstream outstats3;
    std::ifstream testoutstats;
    int filecounter=0;
    std::string tempstr;
    std::ostringstream tempstring;	

    tempstring << filecounter;

    statsfilename += "diffusion_out/diffusion_run";
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

    testoutstats.clear(std::ios::failbit);
    outstats.open(statsfilename.c_str());

    std::cout << statsfilename.c_str() << std::endl;

    std::cout << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " Has_Image: " << image << std::endl;
    std::cout << " # sW = " << s3 << " b = " << s3 - s1 << " mu = " << mu << std::endl;

    outstats << "# Lx: " << lattsize << " Lt: " << nogen << " N_runs: " << numruns << " Has_Image: " << image << std::endl;
    outstats << " # sW = " << s3 << " b = " << s3 - s1 << " mu = " << mu << std::endl;

    cudaStatus = cudaMallocManaged(&lattdown, N * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMallocManaged(&lattup, N * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }
    cudaStatus = cudaMalloc(&devStates, N * sizeof(curandState));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    setup_kernel <<<numBlocks, blockSize>>> (devStates, N, time(0));
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching setup_kernel!\n", cudaStatus);
        goto Error;
    }

    for (int q = 0; q < numruns; q++)
    {
        //initialize lattice
        for (int i = 0; i < N; i++)
        {
            lattdown[i] = 1;
            lattup[i] = 0;
            int j = i / lattsize;
            if (j > left && j <= right) lattdown[i] = 2;
        }

        int t2 = 0;
        int count = 1;
        bool pow_2 = false;
        for (int t = 1; t < nogen; t++)
        {
            if (t == count) {
                pow_2 = true;
                if (t > 1) t2++;
                count = count * 2;
            } else {
                pow_2 = false;
            }

            //even step
            update_odd <<<numBlocks, blockSize>>> (N, lattdown, lattup, lattsize, lattsize2, s1, s3, mu, devStates);
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching update_odd!\n", cudaStatus);
                goto Error;
            }

            //if (pow_2) find_width();

            ++t;
            if (t == count) {
                pow_2 = true;
                if (t > 1) t2++;
                count = count * 2;
            } else {
                pow_2 = false;
            }

            //odd step
            update_even <<<numBlocks, blockSize>>> (N, lattup, lattdown, lattsize, lattsize2, s1, s3, mu, devStates);
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching update_even!\n", cudaStatus);
                goto Error;
            }

            //if (pow_2) find_width();
        }
    }

    int *lattup_host;
    lattup_host = new int[N];
    cudaMemcpy(lattup_host, lattup, N * sizeof(int), cudaMemcpyDeviceToHost);

    if (image)
    {
        for (int i = 0; i < N; i++)
        {
            if (lattup_host[i] == 1) continue;
            int x = i % lattsize;
            int y = i / lattsize;
            double xx = (double)x - 0.5 * (y % 2);
            double yy = (double)y * sqrt(3.0)/2.0;
            outstats << xx << "," << yy << "," << lattup_host[i] << "\n";
        }
        outstats << std::endl;
    }
    
    outstats.close();

    cudaFree(lattdown);
    cudaFree(lattup);

    return 0;

    Error:
        cudaFree(lattdown);
        cudaFree(lattup);
        cudaFree(devStates);
        return (int)cudaStatus;
        
}