#include "Laplacian2DPara.h"

using namespace std ;


//Constructeur :
Laplacian2D::Laplacian2D()
{}
  //Destructeur :
  Laplacian2D::~Laplacian2D()
  {}
    void Laplacian2D::Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double D, double deltaT,\
      int Me, int Np, string Source, string save_all_file)//, string save_points_file, int number_saved_points, vector< vector<double> > &saved_points)
      {
        // // On  initialise les constantes connues de tous les processeurs.
        _x_min = x_min;
        _y_min = y_min;
        _x_max = x_max;
        _y_max = y_max;
        _Nx = Nx;
        _Ny = Ny;
        _D = D;
        _solfinale.resize(_Nx*_Ny);
        _deltaT = deltaT;
        _h_y = (y_max-y_min)/(Ny+1.);
        _h_x = (x_max-x_min)/(Nx+1.);
        _Me = Me;
        _Np = Np;
        _Source = Source;
        _save_all_file = save_all_file;
        // _save_points_file = save_points_file;
        // _number_saved_points = number_saved_points;
        // _saved_points = saved_points;

        system(("rm -Rf "+_save_all_file).c_str());
        system(("mkdir -p ./"+_save_all_file).c_str());

      }

      vector<double> Laplacian2D::InitializeCI(double CI, int i1, int iN, int stencil)
      {
        // // On initialise le vecteur solution ici.
        vector <double> sol_ini(iN-i1+2*stencil);

        for (int i=0; i<iN-i1+2*stencil;i++)
        {
          sol_ini[i] = CI;
        }
        return sol_ini;
      }

      void Laplacian2D::InitializeCL(string gamma_0, string gamma_1)
      {
        // // On initialise les condition limites ici. La configuration quelque peu redondante avec Laplacian2D::Initialize vient d'une ancienne version du code de notre TER dans laquelle on initialisait toutes ces valeurs dans le main et où on ne souhaitait pas avoir trop d'arguments dans la méthode Initialize.
        _gamma_0 = gamma_0;
        _gamma_1 = gamma_1;
      }


      vector<vector<double>> Laplacian2D::InitializeMatrix(int i1, int iN, const int stencil)
      {
        // // On initialise la matrice pentadiagonale ici.
        std::vector<std::vector<double>> LapMat(iN-i1+2*stencil, std::vector<double>(_Nx*_Ny,0.0));
        for (int i=0; i<(iN-i1+2*stencil); i++)
        {
          LapMat[i].resize(_Nx*_Ny);
        }

        double alpha = 1 + 2*_D*_deltaT/(_h_x*_h_x) + 2*_D*_deltaT/(_h_y*_h_y);
        double beta = -_D*_deltaT/(_h_x*_h_x);
        double gamma = -_D*_deltaT/(_h_y*_h_y);

        for (int i=0; i<(iN-i1+2*stencil); i++)
        {
          if (_Me == 0) // i1 = 0 , iN = charge
          {
            LapMat[i][i] = alpha;
            if (i+i1-stencil<_Nx*_Ny-1)
            {
              LapMat[i][i+1] = beta;
            }
            if ( i>0 )
            {
              LapMat[i][i-1] = beta;
            }
            if (i>_Nx-1)
            {
              LapMat[i][i-_Nx] = gamma;
            }
            if (i+_Nx < _Nx*_Ny)
            {
              LapMat[i][i+_Nx] = gamma;
            }
            if (((i+1)%(_Nx) == 0)) // && (i!=0))
            {
              LapMat[i][i+1] = 0.;
            }
            if (((i)%(_Nx) == 0) && (i!=0))
            {
              LapMat[i][i-1] = 0.;
            }
          }

          if (_Me == _Np-1) // i1 = NxNy-1-charge , iN = NxNy-1
          {
            // cout << " je suis la "<< i <<endl;
            // i+=_Nx*_Ny - charge -2*stencil;

            LapMat[i][i+i1-2*stencil] = alpha;
            if (i+i1-2*stencil<_Nx*_Ny-1)
            {
              LapMat[i][i+i1-2*stencil+1] = beta;
            }
            if ( i+i1-2*stencil>0 )
            {
              LapMat[i][i+i1-2*stencil-1] = beta;
            }
            if (i+i1-2*stencil>_Nx-1)
            {
              LapMat[i][i+i1-2*stencil-_Nx] = gamma;
            }
            if (i+i1-2*stencil+_Nx < _Nx*_Ny)
            {
              LapMat[i][i+i1-2*stencil+_Nx] = gamma;
            }

            if (((i +i1 -2*stencil) % _Nx == _Nx-1) && (i!=0))
            {
              LapMat[i][i+i1-2*stencil + 1] = 0.;
            }
            if (((i +i1 -2*stencil ) % _Nx == 0) && (i!=0))
            {
              LapMat[i][i+i1-2*stencil - 1] = 0.;
            }
          }

          if ((_Me != 0) && (_Me != _Np-1))
          {
            LapMat[i][i1-stencil+i]= alpha;
            if (i<_Nx*_Ny-1)
            {
              LapMat[i][i1-stencil+i+1]= beta;
            }
            if ( i+i1-stencil>0 )
            {
              LapMat[i][i1-stencil+i-1]= beta;
            }
            if (i+i1-stencil>_Nx-1)
            {
              LapMat[i][i1-stencil+i-_Nx] = gamma;

            }
            if (i +i1-stencil < _Nx*(_Ny-1))
            {
              LapMat[i][i1-stencil+i+_Nx] = gamma;
            }

            if ((i+i1-stencil+1)%_Nx == 0)
            {
              LapMat[i][i+i1-stencil+1] = 0.;
            }
          }
        }
        return LapMat;
      }

      std::vector<double> Laplacian2D::IterativeSolver(std::vector<std::vector<double>>LapMat, std::vector<double> solloc, int i1, int iN, const int stencil) // le nombre d'itération doit correspondre à la taille du stencil
      {
        // // Cette méthode est au coeur de la résolution du problème, elle permet d'effectuer la marche en temps
        int Nloc = LapMat.size();

        //Sauvegarde d'un points ou plusieurs points particuliers au cours du temps:------------------------------------------------------------------------------------
        ofstream* flux_pts(new ofstream);


        if(_save_points_file != "non")
        {
          //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
          flux_pts->open(_save_points_file+".txt", ios::out);
        }


        //-------------------------------------------------------------------------------------------------------------------------------------------------------------

        for( int i=0 ; i<stencil ; i++)
        {
          //cout << i << " here"<< endl;
          // Laplacian2D::ConditionsLimites(i);
          //--------------------------------------------------------------------------
          //Prise en compte du terme source :
          vector<double> floc(Nloc);
          for (int k=0 ; k < Nloc ;k++)
          {
            floc[k]=0.;
            int num = i1  + k;
            double x = (num%(_Nx))*_h_x;
            double y = (num/(_Nx))*_h_y;

            if (_Source == "non")
            {
              floc[k] = 0+ solloc[k];
            }

            if (_Source == "polynomial")
            {
              floc[k] = 2*_deltaT*(y - y*y + x -x*x) + solloc[k];
            }

            if(_Source == "trigonometrique")
            {
              floc[k] = _deltaT*(sin(x) +cos(y))+ solloc[k];
            }

            if(_Source == "instationnaire")
            {
              floc[k] = _deltaT*exp( -(x/2.)*(x/2.) )*exp( -(y/2)*(y/2) )*cos(M_PI*i*_deltaT/2.)+solloc[k] ;
            }
          }


          //------------------------------------------------------------------------
          for (int k = 0; k < stencil; k++)
          {
            solloc = ConditionsLimites(solloc, i1, iN);
            solloc = CG(LapMat,floc,solloc,0.0001,20, Nloc, Nloc );
          }

        }


        // if ((_save_points_file != "non") and (_Me == 0))
        // {
        //   flux_pts->close();
        // }
        // delete flux_pts;
        //
        // //Barre de chargement
        // if (_Me == 0)//Barre de chargement
        // { printf("\n");}

        return solloc;
      }

      void Laplacian2D::SaveSol(string name_file,std::vector<double> sol)
      {
        // // Cette méthode est celle qui nous permet d'enregistrer la solution sous forme de fichier lisible par paraview. Pour cela, elle envoie tout les vecteurs locaux solloc vers le processeur 0, qui va se charger de reformer le vecteur solution global puis d'écrire le résultat dans le bon format.

        ofstream mon_flux;
        mon_flux.open(name_file, ios::out);
        mon_flux << "# vtk DataFile Version 3.0" << endl
        << "cell" << endl
        << "ASCII" << endl
        << "DATASET STRUCTURED_POINTS" << endl
        << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << endl
        << "ORIGIN 0 0 0" << endl
        << "SPACING " + to_string((_x_max-_x_min)/_Nx)+ " " + to_string((_y_max-_y_min)/_Ny) +" 1" << endl
        << "POINT_DATA " << _Nx*_Ny << endl
        << "SCALARS sample_scalars double" << endl
        << "LOOKUP_TABLE default" << endl;


        for(int i=_Ny-1; i>=0; i--)
        {
          for(int j=0; j<_Nx; j++)
          {
            mon_flux << sol[j + i*_Nx] << " ";
          }
          mon_flux << endl;
        }


        mon_flux.close();
        //Cette partie commentée nous a permis de comparer les résultats théoriques et numériques des solutions pour les cas avec un terme source polynomial et trigonometrique.

      }

      vector<double> Laplacian2D::ConditionsLimites(std::vector<double> solloc, int i1, int iN)
      {
        // // Cette méthode nous permet de mettre à jours le terme source à chaque itération pour prendre en compte les effets des conditions limites.
        // std::cout << "Je suis bien rentre" << std::endl;
        double gamma = -_D*_deltaT/(_h_y*_h_y);
        double beta = -_D*_deltaT/(_h_x*_h_x);

        //////////////////////////CL TEST///////////////////////
        if (_gamma_0 == "non") // HAUT ET BAS
        {
          for (int j = 0; j < _Nx ; j++) //HAUT good
          {
            if ((j<=iN) and (j>=i1))
            solloc[j-i1] = 1;
          }

          for (int j = _Nx*(_Ny-1); j < _Nx*_Ny ; j++) // BAS good
          {
            if ((j<=iN) and (j>=i1))
            solloc[j-i1] = 3;
          }
        }

        if (_gamma_1 == "non") // GAUCHE ET DROITE
        {
          for (int i = 0; i < (_Ny-1)*_Nx; i=i+_Nx) // A GAUCHE good
          {
            if ((i<=iN) and (i>=i1))
            solloc[i -i1] = 2;
          }

          for (int i = _Nx-1; i < _Ny*_Nx; i=i+_Nx) // A DROITE
          {
            if ((i<=iN) and (i>=i1))
            solloc[i -i1] = 4;
          }
        }

        ///////////////////////////////////// CL POLYNOMIAL /////////////////////////////////////
        if (_gamma_0 == "polynomial") // HAUT ET BAS
        {
          for (int j = 0; j < _Nx ; j++) //HAUT good
          {
            if ((j<=iN) and (j>=i1))
            solloc[j-i1] = 0;
          }

          for (int j = _Nx*(_Ny-1); j < _Nx*_Ny ; j++) // BAS good
          {
            if ((j<=iN) and (j>=i1))
            solloc[j-i1] = 0;
          }
        }

        if (_gamma_1 == "polynomial") // GAUCHE ET DROITE
        {
          for (int i = 0; i < (_Ny-1)*_Nx; i=i+_Nx) // A GAUCHE good
          {
            if ((i<=iN) and (i>=i1))
            solloc[i -i1] = 0;
          }

          for (int i = _Nx-1; i < _Ny*_Nx; i=i+_Nx) // A DROITE
          {
            if ((i<=iN) and (i>=i1))
            solloc[i -i1] = 1;
          }
        }

        ///////////////////////////////////// CL TRIGONOMETRIQUE /////////////////////////////////////
        if (_gamma_0 == "trigonometrique") // HAUT ET BAS
        {
          for (int j = 0; j < _Nx ; j++) // En haut
          {
            if ((j<=iN) and (j>=i1))
            {
              solloc[j-i1] = sin(j*_Nx*_h_x) + cos(0);
            }

            for (int j = _Nx*(_Ny-1); j < _Nx*_Ny ; j++) // BAS
            {
              if ((j<=iN) and (j>=i1))
              solloc[j-i1] = (sin(j*_Nx*_h_x) + cos(_Ny*_h_y));
            }
          }

          if (_gamma_1 == "trigonometrique") // GAUCHE ET DROITE
          {
            for (int i = 0; i < (_Ny-1)*_Nx; i=i+_Nx) // A GAUCHE
            {
              if ((i<=iN) and (i>=i1))
              solloc[i -i1] = (sin(0) + cos(i*_h_y));
            }

            for (int i = _Nx-1; i < _Ny*_Nx; i=i+_Nx) // A DROITE
            {
              if ((i<=iN) and (i>=i1))
              solloc[i -i1] = (sin(_h_x*_Nx) + cos(i*_h_y));
            }
          }

          ///////////////////////////////////// CL INSTATIONNAIRE /////////////////////////////////////
          if (_gamma_0 == "instationnaire") // HAUT ET BAS
          {
            for (int j = 0; j < _Nx ; j++) //HAUT
            {
              if ((j<=iN) and (j>=i1))
              solloc[j-i1] = 0;
            }

            for (int j = _Nx*(_Ny-1); j < _Nx*_Ny ; j++) // BAS
            {
              if ((j<=iN) and (j>=i1))
              solloc[j-i1] = 0;
            }
          }

          if (_gamma_1 == "instationnaire") // GAUCHE ET DROITE
          {
            for (int i = 0; i < (_Ny-1)*_Nx; i=i+_Nx) // A GAUCHE
            {
              if ((i<=iN) and (i>=i1))
              solloc[i -i1] = 1;
            }

            for (int i = _Nx-1; i < _Ny*_Nx; i=i+_Nx) // A DROITE
            {
              if ((i<=iN) and (i>=i1))
              solloc[i -i1] = 1;
            }
          }

        }


        return solloc;
      }
