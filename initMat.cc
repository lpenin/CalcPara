// /*
// std::vector<std::vector<double>> Laplacian2D::InitializeMatrix(int i1, int iN, const int stencil)
// {
// */
//
// #include <cmath>
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <algorithm>
// #include <memory>
// #include <vector>
// #include <stdio.h>
// #include <math.h>
//
//
// int main() {
//
//   int _Nx = 3;
//   int _Ny = 3;
//   int _Me = 1; //3; 6;
//   int i1 = 3;//_Me*4;
//   int iN = 5;//(_Me+1)*4;
//   int stencil = 1;
//
//   std::vector<std::vector<double>> LapMat(iN-i1+2*stencil, std::vector<double>(_Nx*_Ny,0.0));
//   double charge = iN-i1;
//   // for (int i=0; i<(charge+2*stencil); i++)
//   // {
//   //   LapMat[i].resize(_Nx*_Ny);
//   // }
//   int Np = 7;
//
//   double alpha = 1;//1 + 2*_D*_deltaT/(_h_x*_h_x) + 2*_D*_deltaT/(_h_y*_h_y);
//   double beta  = 2;//-_D*_deltaT/(_h_x*_h_x);
//   double gamma = 3;//-_D*_deltaT/(_h_y*_h_y);
//
//
//   for (int i=0; i<(charge+2*stencil); i++)
//   {
//     if (_Me == 0) // i1 = 0 , iN = charge
//     {
//       LapMat[i][i] = alpha;
//       if (i+i1-stencil<_Nx*_Ny-1)
//       {
//         LapMat[i][i+1] = beta;
//       }
//       if ( i>0 )
//       {
//         LapMat[i][i-1] = beta;
//       }
//       if (i>_Nx-1)
//       {
//         LapMat[i][i-_Nx] = gamma;
//       }
//       if (i+_Nx < _Nx*_Ny)
//       {
//         LapMat[i][i+_Nx] = gamma;
//       }
//       if (((i+1)%(_Nx) == 0)) // && (i!=0))
//       {
//         LapMat[i][i+1] = 0.;
//       }
//       if (((i)%(_Nx) == 0) && (i!=0))
//       {
//         LapMat[i][i-1] = 0.;
//       }
//
//     }
//
//
//     if (_Me == Np-1) // i1 = NxNy-1-charge , iN = NxNy-1
//     {
//       i+=_Nx*_Ny - charge -2*stencil;
//
//       LapMat[i][i] = alpha;
//       if (i<_Nx*_Ny-1)
//       {
//         LapMat[i][i+1] = beta;
//       }
//       if ( i>0 )
//       {
//         LapMat[i][i-1] = beta;
//       }
//       if (i>_Nx-1)
//       {
//         LapMat[i][i-_Nx] = gamma;
//       }
//       if (i+_Nx < _Nx*_Ny)
//       {
//         LapMat[i][i+_Nx] = gamma;
//       }
//       // if ((i%(_Nx-1) == 0) && (i!=0))
//       // {
//       //   LapMat[i][i+1] = 0.;
//       //   LapMat[i-1][i] = 0.;
//       // }
//     }
//
//
//     if ((_Me != 0) && (_Me != Np-1))
//     {
//       // i+=i1-stencil;
//
//       // LapMat[i][i] = alpha;
//
//       LapMat[i][i1-stencil+i]= alpha;
//
//       if (i<_Nx*_Ny-1)
//       {
//         // LapMat[i][i+1] = beta;
//         LapMat[i][i1-stencil+i+1]= beta;
//       }
//       if ( i+i1-stencil>0 )
//       {
//         // LapMat[i][i-1] = beta;
//         LapMat[i][i1-stencil+i-1]= beta;
//       }
//       if (i+i1-stencil>_Nx-1)
//       {
//         // LapMat[i][i-_Nx] = gamma;
//         LapMat[i][i1-stencil+i-_Nx] = gamma;
//
//       }
//       // if (i+_Nx < _Nx*_Ny)
//       if (i +i1-stencil < _Nx*(_Ny-1))
//       {
//         // LapMat[i][i+_Nx] = gamma;
//         LapMat[i][i1-stencil+i+_Nx] = gamma;
//       }
//       // if (((i+i1-stencil)%(_Nx-1) == 0) && (i!=0))
//       // {
//       //   // LapMat[i][i+1] = 0.;
//       //    // LapMat[i][i+i1-stencil+1] = 0.;
//       //   // LapMat[i-1][i] = 0.;
//       //   // LapMat[i-1][i+i1-stencil] = 0.;
//       // }
//       // if ((((i+i1-stencil)+1)%(_Nx) == 0)) // && (i!=0))
//       // {
//       //   LapMat[i+i1-stencil][i+1] = 0.;
//       // }
//
//       // A changerrrrr en dessous /TODO
//       if (((i+i1-stencil+1)%(_Nx) == 0) && (i+i1-stencil+1<charge+2*stencil))
//       {
//         LapMat[i+i1-stencil+1][i+i1-stencil+2] = 8.;
//         std::cout << i << " " << i << std::endl;
//       }
//     }
//
//
//   }
//   //return LapMat;
//
//   for (int i=0; i<iN-i1+2*stencil;i++ )
//   {
//
//     for (int j=0; j<_Nx*_Ny;j++ )
//     {
//       std::cout << LapMat[i][j] << "   ";
//     }
//     std::cout <<std:: endl;
//   }
//
//   return 0;
// }
