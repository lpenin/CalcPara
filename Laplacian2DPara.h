#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "fonction.h"

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279
#endif


// enum source_enum {non, polynomial, trigonometrique, instationnaire};

class Laplacian2D // pas fini de modifier
{

  protected: // Les attributs de la classe
    double _x_min, _x_max, _y_min, _y_max, _h_x, _h_y, _D, _deltaT;
    int _Nx, _Ny;
    std::vector<std::vector<double> > _LapMatloc; // matrice creuse du laplacien
    std::vector<double> _floc; // vecteur source _f qui prend les données de _sol(i) pour calculer _sol(i+1)
    std::vector<double> _solloc; // vecteur solution U local

    std::string _gamma_0, _gamma_1;

    int _Me,_Np;

    std::string _Source;

    std::string _save_all_file;

    std::string _save_points_file;
    int _number_saved_points;
    std::vector<std::vector <double> > _saved_points;



  public: // Méthodes et opérateurs de la classe

    Laplacian2D();// Constructeur : Initialiser _x_min, _x_max, _y_min; _y_max; _N; _h; _LapMat; _x; _y et _sol.
    virtual ~Laplacian2D();

    std::vector<double> _solfinale;

    void Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double a, double deltaT, int Me, int Np, std::string Source, std::string save_all_file);//, std::string _save_points_file, int number_saved_points,std::vector<std::vector <double> > &saved_points );

    void InitializeCL(std::string gamma_0, std::string gamma_1);

    std::vector<std::vector<double>> InitializeMatrix(int i1,int iN,const int stencil);

    std::vector<double> InitializeCI(double CI, int i1, int iN, int stencil);

    std::vector<double> IterativeSolver(std::vector<std::vector<double>> LapMat, std::vector<double> solloc,int i1, int iN, const int stencil);   // Résout le système _LapMat * _sol = _f avec un solveur itératif.

    void SaveSol(std::string name_file, std::vector<double> sol); // Écrit les solutions dans le fichier "name_file".

    std::vector<double> ConditionsLimites(std::vector<double> solloc, int i1, int iN);
  };



// PRECEDENT INITIALIZATIONS
  //void Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double a, double deltaT, int Me, int Np, string Source, source_enum save_all_file, source_enum _save_points_file, int number_saved_points,std::vector<std::vector <double> > &saved_points );
  //void InitializeCL(source_enum gamma_0, source_enum gamma_1);
