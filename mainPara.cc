#include "Laplacian2DPara.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;


int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);
  MPI_Status status;

  int Me, Np;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);



  std::vector<double> tab_charge(Np);

  int Nx = 4;
  int Ny = 4;
  int i1,iN;

  double xmin = 0.;
  double xmax = 1.;
  double ymin = 0.;
  double ymax = 1.;

  double hx = (xmax-xmin)/(Nx+1);

  double D = 1.;///(1500.*1000.); //Mettre 1. si on fait les cas tests de l'énoncé et 1./(1500.*1000.) si on veut comparer avec notre TER.
  double deltaT =0.9*pow(hx,2)/(4.*D);
  double tfinal = 1;

  std::vector<double> solfinale_loc(Nx*Ny);
  std::vector<double> solfinale(Nx*Ny);

  int stencil=1;
  // if ( Np == 1) {const int stencil=0;}
  // else   {const int stencil=1;}

  double sol_s_deb, sol_s_fin;
  string gamma_0 = "non";
  string gamma_1 = "non";

  int nb_iterations = int(ceil(tfinal/deltaT));

  string Source = "non"; //Peut prendre comMe valeur "non", "polynomial", "trigonoMetrique" ou "instationnaire".
  //Choisir trigonométrique Met à jour les conditions limites automatiqueMent.

  double CI = 0.;


  string save_all_file = "TER"; //Mettre "non" si on ne souhaite pas enregistrer la solution globale au cours du temps sous une forMe lisible par paraview

  charge(Nx*Ny,Np,Me,i1,iN);

  int size_loc= iN-i1+1+2*stencil;

  std::cout << "-----------i1 " << i1 << " iN " << iN << " taille "<< size_loc << std::endl;

  std::vector<double> solloc(size_loc);


  Laplacian2D Lap;
  // Initialisation de toutes les variables, ne pas toucher...
  Lap.Initialize(xmin,xmax,ymin,ymax,Nx,Ny,D,deltaT, Me, Np, Source, save_all_file);//, save_points_file, number_saved_points, saved_points);

  //Initialise _solloc la solution sur un proc + _floc le terMe f sur un proc
  for (int ini=0; ini<Nx*Ny; ini++)
  {
    solfinale_loc[ini]=0.;
    solfinale[ini]=0.;
  }

  solloc=Lap.InitializeCI(CI, size_loc);

  //Sauvegarde des CIs
  for (int j=i1; j<iN+1; j++)
    {
      if (Me == 0)
      {
        solfinale_loc[j]=2; //solloc[j];
      }
      if ((Me != 0) && (Me!=Np-1))
      {
        solfinale_loc[j]=5; //solloc[j-i1+stencil];
      }
      if (Me == Np-1)
      {
        solfinale_loc[j]=10; //solloc[j-i1+2*stencil];
      }
    }

  for (int ini=0; ini <Nx*Ny; ini++)
  {
    MPI_Reduce(&solfinale_loc[ini], &solfinale[ini], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if ((Me==0))
  {
    // for (int i=0; i<Nx*Ny; i++)
    // {
    //   cout<< "    "<<solfinale[i] << endl;
    // }
    Lap.SaveSol("TER/sol_it_"+to_string(0)+".vtk", solfinale); // -> partage solloc avec le processeur 0 pour qu'il puisse écrire la solution dans un fichier .vtk
  }
  //Fin sauvegarde

  // Initialise les CLs
  Lap.InitializeCL(gamma_0, gamma_1);

  // Initialise la petite Matrice par proc
  vector<vector<double>> LapMat(size_loc);

  for (int i=0; i<(size_loc); i++)
  {
    LapMat[i].resize(Nx*Ny);
  }

  LapMat=Lap.InitializeMatrix(i1, iN, stencil, size_loc);

  auto start = chrono::high_resolution_clock::now();

  // Nouvelle boucle sur le nombre d'itétration
  for (int i=0; i<2; i+= stencil) //TODO
  // for (int i=0; i<nb_iterations; i+= stencil)
  {
    // if (Me==0)
    // {
    //   for (int test=0; test<iN-i1+2; test++)
    //   {
    //     cout <<"Me= " << Me <<" pour i= "<< i1+test << " solloc = "<< solloc[test]<< endl;
    //   }
    // }

    solloc=Lap.IterativeSolver(LapMat, solloc, i1, iN, stencil);
    // cout << " ICI encore me = " << Me << " i1, iN "<< i1 << " et "<<iN << " solloc size= "<< solloc.size()<<endl;
    if ((Me==0))
    {
      for (int i=0; i<size_loc; i++)
      {
        cout<<"  "<< Me <<"   "<<solloc[i] << endl;
      }
    }

    if ((Me==1))
    {
      for (int i=0; i<size_loc; i++)
      {
        cout<<"  "<< Me <<"   "<<solloc[i] << endl;
    }
  }

    if ((Me==2))
    {
      for (int i=0; i<size_loc; i++)
      {
        cout<<"  "<< Me <<"   "<<solloc[i] << endl;
    }
  }

    if ((Me==3))
    {
      for (int i=0; i<size_loc; i++)
      {
        cout<<"  "<< Me <<"   "<<solloc[i] << endl;
    }
  }


    //communication condition proc stencil
    for (int ini=0; ini <Nx*Ny; ini++)
    {
      solfinale_loc[ini]=0.;
      solfinale[ini]=0.;
    }

    for (int k=1; k< Np-1; k++)
    {
      for (int s=0; s<stencil; s++)
      {
        // if ((Me == k) && (Me != 0) && (Me!=Np))
        if (Me == k)
        {
          if (Me== 1){MPI_Send(&solloc[2*stencil+s], 1, MPI_DOUBLE, k-1, 1997, MPI_COMM_WORLD);}
          else {MPI_Send(&solloc[stencil+s], 1, MPI_DOUBLE, k-1, 1997, MPI_COMM_WORLD);}

          if (Me== Np-1){MPI_Send(&solloc[solloc.size()-1-3*stencil+s], 1, MPI_DOUBLE, k+1, 1997, MPI_COMM_WORLD);}
          else{MPI_Send(&solloc[solloc.size()-1-2*stencil+s], 1, MPI_DOUBLE, k+1, 1997, MPI_COMM_WORLD);}


          MPI_Recv(&sol_s_deb, 1, MPI_DOUBLE, k-1, 1997, MPI_COMM_WORLD, &status);
          solloc[s]=sol_s_deb;
          MPI_Recv(&sol_s_fin, 1, MPI_DOUBLE, k+1, 1997, MPI_COMM_WORLD, &status);
          solloc[solloc.size()-stencil+s]=sol_s_fin;

          for (int j=i1; j<iN+1;j++)
          {
            solfinale_loc[j]=solloc[j-i1+stencil];
          }
        }
      }
    }

    for (int s=0; s<stencil; s++)
    {
      if (Me ==0)
      {
        MPI_Send(&solloc[solloc.size()-3*stencil+s], 1, MPI_DOUBLE, 1, 1997, MPI_COMM_WORLD);
        MPI_Recv(&sol_s_fin, 1, MPI_DOUBLE, 1, 1997, MPI_COMM_WORLD, &status);
        solloc[solloc.size()-stencil+s]=sol_s_fin;
        for (int j=i1; j<iN+1;j++)
        {
          solfinale_loc[j]=solloc[j-i1];
        }
      }
      if (Me == Np-1)
      {
        MPI_Send(&solloc[2*stencil+s], 1, MPI_DOUBLE, Np-2, 1997, MPI_COMM_WORLD);
        MPI_Recv(&sol_s_deb, 1, MPI_DOUBLE, Np-2, 1997, MPI_COMM_WORLD, &status);
        solloc[s]=sol_s_deb;
        for (int j=i1; j<iN+1;j++)
        {
          solfinale_loc[j]=solloc[j-i1+2*stencil];
        }
      }
    }

    //Barrière de communications
    for (int ini=0; ini <Nx*Ny+1; ini++)
    {
      MPI_Reduce(&solfinale_loc[ini], &solfinale[ini], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD); // Useless post reduction ?


    //Barre de chargeMent
    if (Me == 0)
    {
      int i_barre;
      int p = floor((((double)i)/((double)nb_iterations))*100);
      printf( "[" );
      for(i_barre=0;i_barre<=p;i_barre+=2) printf( "*" );
      for (;i_barre<=100; i_barre+=2 ) printf( "-" );
      printf( "] %3d %%", p );

      for(i_barre=0;i_barre<59;++i_barre) printf( "%c", 8 );

      fflush(stdout );
    }

    //  if (Me==0)
    if ((Me==0) and (i%10==0))
    {
      int iter_num = i/10 + 1;
      Lap.SaveSol("TER/sol_it_"+to_string(iter_num)+".vtk", solfinale); // -> partage solloc avec le processeur 0 pour qu'il puisse écrire la solution dans un fichier .vtk
    }

  }


  // Lap.IterativeSolver(nb_iterations);
  auto finish = chrono::high_resolution_clock::now();

  double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();

  if(Me ==0)
  {
    cout << endl << "Le prog a mis " << t*0.000001 << " secondes a s'effectuer" << endl;
  }
  MPI_Finalize();

  return 0;
}
