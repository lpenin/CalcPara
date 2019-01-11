#include "fonction.h"


void charge(int N, int Np, int Me, int& i1, int& iN)
{
  int q = N/Np;
  int r = N - q*Np;

  if (Me < r)
  {
    i1 = (q+1)*Me;
    iN = (q+1)*(Me+1)-1;
  }
  else
  {
    i1 = q*Me + r;
    iN = q*(Me+1) + r - 1;
  }
}

// Matrix times vector
std::vector<double> prodMatVec( std::vector<std::vector<double>>A, std::vector<double> V )
{
  int n = A.size();
  std::vector<double> C( n );
  for ( int i = 0; i < n; i++ ) C[i] = dot( A[i], V );
  return C;
}

double dot(std::vector<double> U, std::vector<double> V)
{
  int N = U.size();
  double y=0;
  for (int i=0; i<N; i++)
  {
    y += U[i]*V[i];
  }
  return y;
}




// Gradient Conjugué
std::vector<double> CG(std::vector<std::vector<double> > Aloc, std::vector<double> bloc, std::vector<double> x0loc , double err, int kmax, int nx, int ny)
{
  int k;

  double norm_r, nr_carre, nr2_carre;
  std::vector<double> wloc, rloc, r_1loc, ploc, dloc, xloc;

  int Np = nx;

  wloc.resize(Np); rloc.resize(Np); ploc.resize(Np); dloc.resize(Np);

  xloc = x0loc;


  wloc = prodMatVec(Aloc, xloc);

  for (int i = 0; i < Np; i++)
  {
    rloc[i] = bloc[i] - wloc[i];
    ploc[i] = rloc[i];
    // cout << rloc[i] << "  rloc GROOOOSS"<< endl;
  }

  k = 0;

  nr_carre = dot(rloc,rloc);
  norm_r = sqrt(nr_carre);
  // cout <<" norm r "<< norm_r << " et eps "<< err<<endl;
  while ((norm_r > err) and (k<kmax))
  {

    dloc = prodMatVec(Aloc, ploc);
    double alpha = nr_carre/dot(ploc,dloc);

    for (int i = 0; i < Np; i++)
    {
      xloc[i] += alpha*ploc[i];
      rloc[i] -= alpha*dloc[i]; //rk devient r(k+1)
    }

    nr2_carre = dot(rloc,rloc); //norme de r(k+1)

    double beta = nr2_carre/nr_carre;

    for (int i = 0; i < Np; i++)
    {
      ploc[i] = rloc[i] + beta*ploc[i];
    }

    nr_carre = nr2_carre;
    norm_r = sqrt(nr_carre);

    k++;
  }

  return xloc;

}

void printvect(std::vector<double> uloc)
{
  // // Fonction qui permet d'afficher les composante d'un vecteur global en gardant la configuration locale. Cette fonction n'est plus utilisée dans le code mais nous a permis de le débugger.


  int Me, Np;
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  for (int i =0; i < Np; i++)
  {
    if (Me == i)
    {
      std::cout << "Me = " << i << std::endl;
      for (int j = 0; j < uloc.size(); j++)
      {
        std::cout << uloc[j] << std::endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (Me == Np-1)
  {
    std::cout << "fin affichage" << std::endl;
  }
}
