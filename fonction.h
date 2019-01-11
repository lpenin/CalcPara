#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>

void charge (int N, int Np, int Me, int& i1, int& iN);

std::vector<double> prodMatVec(std::vector<std::vector<double>>A, std::vector<double> V);

double dot(std::vector<double> U, std::vector<double> V);

std::vector<double> CG (std::vector<std::vector<double> > Aloc, std::vector<double> bloc, std::vector<double> x0loc, double err, int kmax, int nx, int ny);
void printvect( std::vector<double> uloc);
