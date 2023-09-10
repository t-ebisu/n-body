#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <fstream>
#include <strstream>
#include <sstream> //io
#include "gettime.h"
using namespace std;


static double Mass = 1.0;
static double G = 1.0;
static double Radius = 1.0;
static double eps = 0.03125;
static int max_n = 16384;
static int max_nstep = 100000;


typedef struct Particle{
  double pos[3];
  double vel[3];
  double acc[3];
  double pot;
  double mass;
  //long long id;
}Particle, *pParticle;



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


double calc_vrms(double r_vir, pParticle particle, int n){
  double w = 0.0;
  for(int i=0; i<n; i++) w += particle[i].pot; 
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
double calc_E(pParticle particle, int n){
  double k = 0.0, w = 0.0;
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++) k += particle[i].mass * particle[i].vel[j] * particle[i].vel[j];
    w -= particle[i].pot;
  }

  cout << "k:" << k/2 << ", w:" << w/2 << endl;

  return 0.5*(k+w);
}


//Calculate gravity
void calc_grav(pParticle particle, double eps2, int n){
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++) particle[i].acc[j] = 0.0;
    particle[i].pot = 0.0;
  }
  
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      double r[3];
      for(int k=0; k<3; k++) r[k] = (particle[j].pos[k] - particle[i].pos[k]);
      
      double r2 = eps2;
      for(int k=0; k<3; k++) r2 += r[k]*r[k];
      
      double rinv = 1.0 / sqrt(r2);
      double r3inv = rinv * rinv * rinv;
      
      for(int k=0; k<3; k++){ 
        particle[i].acc[k] += G*particle[j].mass*r[k]*r3inv;
        particle[j].acc[k] -= G*particle[i].mass*r[k]*r3inv;
      }

      particle[i].pot += G*particle[i].mass*particle[j].mass*rinv;
      particle[j].pot += G*particle[i].mass*particle[j].mass*rinv;
    }
  }
}


void runge_kutta4(pParticle particle, double eps2, double dt, int n){
  double xtmp[n][3];
  double k1[n][3], k2[n][3], k3[n][3], k4[n][3];
  double l1[n][3], l2[n][3], l3[n][3], l4[n][3];

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k1[i][j] = dt * particle[i].vel[j];
      l1[i][j] = dt * particle[i].acc[j];
      xtmp[i][j] = particle[i].pos[j];
      particle[i].pos[j] = xtmp[i][j] + k1[i][j]*0.5;
    }
  }

  calc_grav(particle, eps2, n);
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k2[i][j] = dt * (particle[i].vel[j] + l1[i][j]*0.5);
      l2[i][j] = dt * particle[i].acc[j];
      particle[i].pos[j] = xtmp[i][j] + k2[i][j]*0.5;
    }
  }

  calc_grav(particle, eps2, n);
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k3[i][j] = dt * (particle[i].vel[j] + l2[i][j]*0.5);
      l3[i][j] = dt * particle[i].acc[j];
      particle[i].pos[j] = xtmp[i][j] + k3[i][j];
    }
  }

  calc_grav(particle, eps2, n);
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      k4[i][j] = dt * (particle[i].vel[j] + l3[i][j]);
      l4[i][j] = dt * particle[i].acc[j];
      particle[i].pos[j] = xtmp[i][j];
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      particle[i].pos[j] += (k1[i][j] + 2.0*k2[i][j] + 2.0*k3[i][j] + k4[i][j])/6.0;
      particle[i].vel[j] += (l1[i][j] + 2.0*l2[i][j] + 2.0*l3[i][j] + l4[i][j])/6.0;
    }
  }

  calc_grav(particle, eps2, n);
}


void leap_frog(pParticle particle, double eps2, double dt, int n){
  double v_h[n][3];

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      v_h[i][j] = particle[i].vel[j] + 0.5 * particle[i].acc[j] * dt;
      particle[i].pos[j] = particle[i].pos[j] + v_h[i][j] * dt;
    }
  }

  calc_grav(particle, eps2, n);

  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++){
      particle[i].vel[j] = v_h[i][j] + 0.5 * particle[i].acc[j] * dt;
    }
  }
}


//begin io functions
void dump_file(pParticle particle, int n,int nstep, double error){
  string filename = "snapshot_" + to_string(nstep) + "-nstep.log";
  std::cout << "dump " << filename << endl;
  
  ofstream ofs(filename);
  ofs << n << endl;
  ofs << nstep << endl;
  ofs << error << endl;
  for(int i=0; i<n; i++){
    ofs << setprecision(15)
     << setiosflags(ios::scientific)
     << particle[i].pos[0] << "," << particle[i].pos[1] << "," << particle[i].pos[2] << ","
     << particle[i].vel[0] << "," << particle[i].vel[1] << "," << particle[i].vel[2] << ","
	   << particle[i].mass << endl;
  }

  ofs.close();
}


void read_file(string filename, pParticle particle, int *n, int *nstep, double *error){
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

    for(int i=0; i<3; i++){
      getline(stream, data, ',');
      particle[itr].pos[i] = stod(data);
    }
    
    for(int i=0; i<3; i++){
      getline(stream, data, ',');
      particle[itr].vel[i] = stod(data);
    }

    getline(stream, data, ',');
    particle[itr].mass = stod(data);

    itr++;
    // if over max_n or n ; break;
    if(itr >= max_n) break;

  }
  
}
// end io functions


void init(pParticle particle, double r_vir, double eps2, const int n){
  for(int i=0; i<n; i++){
    particle[i].mass = Mass/n;
    pos(particle[i].pos);
  }

  calc_grav(particle, eps2, n);
  double vrms = calc_vrms(r_vir, particle, n);
  cout << "vrms:" << vrms << endl;
  for(int i=0; i<n; i++){
    for(int j=0; j<3; j++) particle[i].vel[j] = vrms * gaussian();
  }
}



int  main(){

  int n;
  int nstep = 0;
  double tnow = 0, dt, tend;
  double eps2 = eps * eps;
  pParticle particle = new Particle[max_n];
  double e_0, e_end;
  double r_vir;
  
  //set default params
  n = 128; //number of particles
  dt = pow(2,-5);
  tend = 1.0;
  r_vir = 0.5; //virial radius
  //string filename = "hogehoge.csv";

  //cin << filename;
  //or
  //cin << n << dt << tend << r_vir; 

  srand48(1);
  init(particle, r_vir, eps2, n);
  //read_file(filename, particle, &n, &nstep, &e_0);
  
  double nowtime = 0.0;
  getTime(&nowtime);

  e_0 = calc_E(particle, n);
  
  while(tnow < tend){
    //select runge_kutta or leap_frog
    runge_kutta4(particle, eps2, dt, n);
    //leap_frog(particle, eps2, dt, n);

    tnow += dt;
    nstep++;
    
    if(nstep % 10 == 0){
      e_end = calc_E(particle, n);
      dump_file(particle, n, nstep, fabs((e_end-e_0)/e_0));
    }

  }

  double difftime = getTime(&nowtime);

  e_end = calc_E(particle, n);
  cout << scientific << "Error : " << fabs((e_end-e_0)/e_0) << endl;
  cout << scientific << "time : " << difftime << endl;

  delete[] particle;

  return 0;
}
