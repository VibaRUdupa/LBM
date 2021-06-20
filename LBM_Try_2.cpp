// conditions and one wall moving with tangential velocity of 6 m/s.
// Viba R Udupa, June 19, 2021

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define nx 101
#define ny 101

//// Parameter variables initialization
float SL, LU, kv, Re; // Square edge length, Lid speed, kinematic viscosity, Reynolds number
float l_LU, l_kv, dx, dy, tau; // Lattice - lid speed, kinematic viscosity, distance between nodes and relaxation time.

int r, c, a, a1, ra, ca;    // Looping variables i for x, j for y and a for the 9 vectors.
int is_wall[nx][ny];  //Indicator variable of if the node is a wall node (=1)

float f[9][nx][ny], ft[9][nx][ny], feq[9][nx][ny]; // f stores currect distribution and f_temp is a temporary distribution storage after collision step.
float wt[9] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; //Stores the weigth of every direction associated with D2Q9
float ux, uy, ro, rho_avg, u2, term1, term2; //Stores the values of velocity and density in lattice units.
float rho[nx][ny];
float u[nx][ny], v[nx][ny];

int ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
int ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
// ex and ey stores the coordinates of the directions.
int time = 20000, ts; // Run for 20,000 time steps.

float rho_0; // Initial density is taken as 1 units
float sum, usum, vsum;

FILE *f1, *f2, *f3, *f4, *f5;

// User defined Functions
void problem_spec();
void initialization();
void time_step(int ts);
void distribution_collision();
void streaming();
void Boundary_Conditions();
void Calculations();
void Results();


int main()
{
  problem_spec();
  initialization();
  printf("\nProblem specified\nInitialization complete\nIn main\n");
  for(ts = 1; ts<=time; ts++)
  {
    time_step(ts);
  }

}

void problem_spec()
{
  SL = 0.2;
  LU = 6.0;
  kv =  0.0012;
  Re = LU*SL/kv;

  l_LU = 0.1;
  l_kv = 0.01;
  dx = 1;
  dy = 1;
  tau = 3*l_kv + 0.5;
  rho_0 = 1.0;
  //printf("\nLBM for Lid Driven Cavity\nl_LU = %f\nl_kv = %f\nRe = %f\ntau = %f\n", l_LU, l_kv, l_LU*(nx-1)/l_kv, tau);
}

void initialization()
{
  // Initialize the initial distribution function
  for ( c = 0; c < ny; c++)
  {
    for (r = 0; r < nx; r++)
    {
      rho[r][c] = rho_0;
      u[r][c] = 0.0;
      v[r][c] = 0.0;
    }
  }

  for(c = 1; c < ny-1; c++)
  {
    u[0][c] = l_LU;
    v[0][c] = 0.0;
  }

  //File openings
  f1 = fopen("lcd_uvfield.txt", "w");
  f2 = fopen("lcd_uvely.txt", "w");
  f3 = fopen("lcd_vvelx.txt", "w");
  f4 = fopen("lcd_timeu.txt", "w");
  f5 = fopen("lcd_streamf.txt", "w");
}

void time_step(int ts)
{
  distribution_collision();
  streaming();
  Boundary_Conditions();
  Calculations();
  if(ts % 200 == 0)
    printf("\nt = %d\t %f\t%f\t%f\n", ts, u[(nx-1)/2 +1][(ny-1)/2 +1],v[(nx-1)/2 +1][(ny-1)/2 +1], rho[(nx-1)/2 +1][(ny-1)/2 +1]);
  //Results();
}

void distribution_collision()
{
  for (r = 0; r < nx ; r++)
  {
    for (c = 0; c < ny ; c++)
    {
      term1 = u[r][c]*u[r][c] + v[r][c]*v[r][c];
      for ( a = 0 ; a < 9 ; a ++ )
      {
        term2 = u[r][c]*ex[a] + v[r][c]*ey[a];
        feq[a][r][c] = wt[a]*rho[r][c]*(1.0 + 3.0*term2 + 4.5*term2*term2 - 1.5*term1);
        f[a][r][c] = f[a][r][c] - (f[a][r][c] - feq[a][r][c])/tau;
      }
    }
  }
  //printf("\nt = %4d \tf[4][50][50] =%f \tf[4][51][51] = %f",ts, f[4][50][50],f[4][51][51] );
}

void streaming()
{
  for(r = 0; r <= ny-1 ; r++)
  {
    for(c = nx-1; c >=1 ; c--)
    {
      a = 1;
      f[a][r][c] = f[a][r][c-1];
    }

    for(c = 0; c < nx-1 ; c++)
    {
      a = 3;
      f[a][r][c] = f[a][r][c+1];
    }
  }

  for(r = nx-1; r >= 1 ; r--)
  {
    for(c = 0; c <= nx-1 ; c++)
    {
      a = 2;
      f[a][r][c] = f[a][r-1][c];
    }

    for(c = nx-1; c >= 1 ; c--)
    {
      a = 5;
      f[a][r][c] = f[a][r-1][c-1];
    }

    for(c = 0; c < nx-1 ; c++)
    {
      a = 6;
      f[a][r][c] = f[a][r-1][c+1];
    }
  }

  for(r = 0; r <ny-1 ; r++)
  {
    for(c = 0; c <= nx-1 ; c++)
    {
      a = 4;
      f[a][r][c] = f[a][r+1][c];
    }

    for(c = nx-1; c >= 1 ; c--)
    {
      a = 8;
      f[a][r][c] = f[a][r+1][c-1];
    }

    for(c = 0; c < nx-1 ; c++)
    {
      a = 7;
      f[a][r][c] = f[a][r+1][c+1];
    }
  }
}

void Boundary_Conditions()
{
  for(r=0 ; r < ny; r++)
  {
      // West boundary bounceback
      f[1][r][0] = f[3][r][0];
      f[5][r][0] = f[7][r][0];
      f[8][r][0] = f[6][r][0];

      // East boundary bounceback
      f[3][r][nx-1] = f[1][r][nx-1];
      f[7][r][nx-1] = f[5][r][nx-1];
      f[6][r][nx-1] = f[8][r][nx-1];
    }

    for(c=0 ; c < nx; c++)
    {
      // Bottom boundary bounceback
      f[2][ny-1][c] = f[4][ny-1][c];
      f[5][ny-1][c] = f[7][ny-1][c];
      f[6][ny-1][c] = f[8][ny-1][c];

      if(c == 0 || c == ny-1)
      continue;

      // Top moving wall Boundary
      float rh = f[0][0][c] + f[1][0][c] + f[3][0][c] + 2 * (f[2][0][c] + f[6][0][c] + f[5][0][c]);
      f[4][0][c] = f[2][0][c];
      f[8][0][c] = f[6][0][c] + rh*l_LU/LU;
      f[7][0][c] = f[5][0][c] - rh*l_LU/LU;
    }
}

void Calculations()
{
  for(c = 0; c < ny; c++)
  {
    for (r = 0; r< nx; r++)
    {
      sum = 0.0;
      for (a = 0; a<9; a++)
      {
        sum = sum + f[a][r][c];
      }
      rho[r][c] = sum;
    }
  }

  for(c = 1; c < ny; c++)
  {
      rho[0][c] = f[0][0][c] + f[1][0][c] + f[3][0][c] + 2 * (f[2][0][c] + f[6][0][c] + f[5][0][c]);
  }

  for(c = 1; c < ny; c++)
  {
    for(r = 1; r < nx-1; r++)
    {
      usum = 0.0;
      vsum = 0.0;

      for(a=0; a<9; a++)
      {
        usum = usum + f[a][r][c]*ex[a];
        vsum = vsum + f[a][r][c]*ey[a];
      }
      u[r][c] = usum/rho[r][c];
      v[r][c] = vsum/rho[r][c];
    }
  }
  //printf("\trho[51][51] =%f \tu[51][51] = %f \tv[51][51] = %f", rho[51][51],u[51][51],v[51][51]);
}

void Results()
{
  //for(c = 1; c < nx; c++)
  {
    //rho_avg = 0.5 * (rho[nx-1][c-1] + rho[nx-1]][c]);
  }

}
