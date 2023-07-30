#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <fstream>
#include <strstream>
#include "gettime.h"
using namespace std;


static double Mass = 1.0;
static double G = 1.0;
static double Radius = 1.0;
static double eps = 0.03125;
static int max_n = 16384;
static int max_nstep = 100000;


// Gaussian with mean = 0.0 and dispersion = 1.0 by Box-Muller method 
double gaussian(void){
  double x, y, z, r2;
  do{
    x = 2.0*drand48() - 1.0;
    y = 2.0*drand48() - 1.0;
    r2 = x*x + y*y;
  }while(r2 >= 1.0 || r2 == 0.0);
  z = sqrt(-2.0*log(r2)/r2)*x;
  return (z);
}


//Calculate velocity dispersion
double calc_vrms(double r_vir, double *p, int n){
  double w = 0.0;
  for(int i=0; i<n; i++) w += p[i]; 
  return sqrt((r_vir*w)/(3*Mass));
}


void pos(double x[3]){
  double r2;
  do{
    r2 = 0.0;
    for(int i=0; i<3; i++){
      x[i] = Radius * (2.0*drand48() - 1.0);
      r2 += x[i]*x[i];
    }
  }while(r2 > Radius);
}


//Calculate total enargy
double calc_E(double (*x)[3], double (*v)[3], double *m, double *p, int n){
  double k = 0.0, w = 0.0;
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++) k += m[i] * v[i][j] * v[i][j];
    w -= p[i];
  }

  cout << "k:" << k/2 << ", w:" << w/2 << endl;

  return 0.5*(k+w);
}


//Calculate gravity
void calc_grav(double (*x)[3], double *m, double (*a)[3], double *p, double eps2, int n){
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++) a[i][j] = 0.0;
    p[i] = 0.0;
  }
  
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      double r[3];
      for(int k=0; k<3; k++) r[k] = (x[j][k] - x[i][k]);
      
      double r2 = eps2;
      for(int k=0; k<3; k++) r2 += r[k]*r[k];
      
      double rinv = 1.0 / sqrt(r2);
      double r3inv = rinv * rinv * rinv;
      
      for(int k=0; k<3; k++){ 
        a[i][k] += G*m[j]*r[k]*r3inv;
        a[j][k] -= G*m[i]*r[k]*r3inv;
      }

      p[i] += G*m[i]*m[j]*rinv;
      p[j] += G*m[i]*m[j]*rinv;
    }
  }
}


void runge_kutta4(double (*x)[3], double (*a)[3], double (*v)[3], double *m, double *p, double eps2, double dt, int n){
  double xtmp[n][3];
  double k1[n][3], k2[n][3], k3[n][3], k4[n][3];
  double l1[n][3], l2[n][3], l3[n][3], l4[n][3];

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k1[i][j] = dt * v[i][j];
      l1[i][j] = dt * a[i][j];
      xtmp[i][j] = x[i][j] + k1[i][j]*0.5;
    }
  }

  calc_grav(xtmp, m, a, p, eps2, n);
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k2[i][j] = dt * (v[i][j] + l1[i][j]*0.5);
      l2[i][j] = dt * a[i][j];
      xtmp[i][j] = x[i][j] + k2[i][j]*0.5;
    }
  }

  calc_grav(xtmp, m, a, p, eps2, n);
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k3[i][j] = dt * (v[i][j] + l2[i][j]*0.5);
      l3[i][j] = dt * a[i][j];
      xtmp[i][j] = x[i][j] + k3[i][j];
    }
  }

  calc_grav(xtmp, m, a, p, eps2, n);
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k4[i][j] = dt * (v[i][j] + l3[i][j]);
      l4[i][j] = dt * a[i][j];
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      x[i][j] += (k1[i][j] + 2.0*k2[i][j] + 2.0*k3[i][j] + k4[i][j])/6.0;
      v[i][j] += (l1[i][j] + 2.0*l2[i][j] + 2.0*l3[i][j] + l4[i][j])/6.0;
    }
  }

  calc_grav(x, m, a, p, eps2, n);
}


void leap_frog(double (*x)[3], double (*a)[3], double (*v)[3], double *m, double *p, double eps2, double dt, int n){
  double v_h[n][3];

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      v_h[i][j] = v[i][j] + 0.5 * a[i][j] * dt;
      x[i][j] = x[i][j] + v_h[i][j] * dt;
    }
  }

  calc_grav(x, m, a, p, eps2, n);

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      v[i][j] = v_h[i][j] + 0.5 * a[i][j] * dt;
    }
  }
}


//begin io functions
void dump_file(double (*x)[3], double (*v)[3], const int n, const int nstep, double error){
  string fname = "snapshot_" + to_string(nstep) + "-nstep.log";
  cout << "dump " << fname << endl;
  
  ofstream ofs(fname);
  ofs << n << endl;
  ofs << nstep << endl;
  ofs << error << endl;
  for(int i=0; i<n; i++){
    ofs << setprecision(15)
     << setiosflags(ios::scientific)
     << x[i][0] << "," << x[i][1] << "," << x[i][2] << ","
     << v[i][0] << "," << v[i][1] << "," << v[i][2] << ","
	   << endl;
  }

  ofs.close();
}


void read_file(string filename, double (*x)[3], double (*v)[3], double (*m), int *n, int *nstep, double *error){
  ifstream ifs(filename);
  string line, data;

  getline(ifs, line);
  *n = stoi(line);

  getline(ifs, line);
  *nstep = stoi(line);
  
  getline(ifs, line);
  *error = stod(line);
  
  int itr = 0;
  while(getline(ifs, line)){  
    istringstream stream(line);
    /*while(getline(stream, data, ',')){
      //hogehoge
    }*/

    for(int i=0; i<3; i++){
      getline(stream, data, ',');
      x[itr][i] = stod(data);
    }
    
    for(int i=0; i<3; i++){
      getline(stream, data, ',');
      v[itr][i] = stod(data);
    }

    getline(stream, data, ',');
    m[itr] = stod(data);

    itr++;
    // if over max_n or n ; break;

  }
  
}
// end io functions


//Initialize information of the particles
void init( double (*x)[3], double (*v)[3], double *m, double (*a)[3], double *p, double r_vir, double eps2, const int n){
  for(int i=0; i<n; i++){
    m[i] = Mass/n;
    pos(x[i]);
  }

  calc_grav(x, m, a, p, eps2, n);
  double vrms = calc_vrms(r_vir, p, n);
  cout << "vrms:" << vrms << endl;
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++) v[i][j] = vrms * gaussian();
  }
}



int  main(){

  int n;
  int nstep = 0;
  double tnow = 0, dt, tend;
  double eps2 = eps * eps;
  double x[max_n][3];
  double v[max_n][3];
  double a[max_n][3];
  double p[max_n];
  double m[max_n];
  double e_0, e_end;
  double r_vir;

  n = 20; //number of particles
  dt = pow(2,-5);
  tend = 1.0;
  r_vir = 0.5; //virial radius
  

  srand48(1);
  init(x, v, m, a, p, r_vir, eps2, n);
  //read_file(filename, x, v, m, &n, &nstep, &e_0);

  double nowtime = 0.0;
  getTime(&nowtime);

  e_0 = calc_E(x, v, m, p, n);
  
  while(tnow < tend){
    //runge_kutta4(x, a, v, m, p, eps2, dt, n);
    leap_frog(x, a, v, m, p, eps2, dt, n);

    tnow += dt;
    nstep++;
    
    if(nstep % 10 == 0){
      e_end = calc_E(x, v, m, p, n);
      dump_file(x, v, n, nstep, fabs((e_end-e_0)/e_0));
    }

    //for debug
    //fprintf(stderr, "%e, %e, %e, %e, %e, %e, %e, %e, %e\n", x[0][0], x[0][1], x[0][2], v[0][0], v[0][1], v[0][2], a[0][0], a[0][1], a[0][2]);

  }

  double difftime = getTime(&nowtime);

  e_end = calc_E(x, v, m, p, n);
  cout << scientific << "Error : " << fabs((e_end-e_0)/e_0) << endl;
  cout << scientific << "time : " << difftime << endl;

  return 0;
}
