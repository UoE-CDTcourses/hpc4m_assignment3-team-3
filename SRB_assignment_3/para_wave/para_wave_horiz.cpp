#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

static const double pi = 3.1415926536;
static const int M = 2306;
static const double T = 1;

int main(int argc, char* argv[]){

   int rank, nproc, ierr, rank1, rank2;
   double dt = 0.2/(M);
   int N = T/dt;
   int t1=0.333/dt, t2=0.666/dt, t3=N; // points at which we want to print our results.
   double dy = 2./(M-1);
   double dx = dy;
   double dy2 = dy*dy;
   double dtdy2 = dt*dt/dy2;
   double x,y;
   MPI_Comm comm;

   MPI_Init(NULL,NULL);

   MPI_Barrier(MPI_COMM_WORLD);
   double t_start = MPI_Wtime();

   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   int J = 2+(M-2)/nproc;
   double to_send[M-2],to_recv[M-2];

   double U[3][J][M];
   double U_final[M][M];
   double U_recv[J][M];

 
 
   // Functional initial condition.
       for (int i=0; i<J; ++i){
         for (int j=0; j<M; ++j){
            y = (rank*(J-2)+i)*dy - 1;
            x = j*dx-1;

            U[0][i][j] = U[1][i][j] = exp( -40 *  ((x-0.4)*(x-0.4) + (y*y )));
         }
      } 
   // left and right boundaries needed for every rank.
    for(int i=0;i<J;i++){
      U[0][i][0] = U[0][i][J-1] = 0;
      U[1][i][0] = U[1][i][J-1] = 0;
      U[2][i][0] = U[2][i][J-1] = 0;
   } 
   // Top and bottom boundaries required for the first and last rank.

    if(rank == 0){
      for(int j = 1;j<M-1;j++){
         U[0][0][j] = U[1][0][j]= U[2][0][j] = 0;
      }
   }
   if(rank == nproc-1){
      for(int j = 1;j<M-1;j++){
         U[0][J-1][j] = U[1][J-1][j]= U[2][J-1][j] = 0;
      }
   } 
   // Printing initial conditions 
  if(rank !=0){
     MPI_Ssend(U,M*J,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
  }

  if(rank == 0){
     for(int i = 0;i<J-2;i++){
        for(int j =0; j<M;j++){
           U_final[i][j] = U[0][i][j];
        }
     }
     for(int send_rank = 1;send_rank<nproc;send_rank++){
        MPI_Recv(U_recv,M*J,MPI_DOUBLE,send_rank,send_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        for(int i =0; i<J-2;i++){
           for(int j=0;j<M;j++){
            U_final[send_rank*(J-2)+i][j] = U_recv[i][j];
           }
        }
        if(send_rank == nproc-1){
           for(int i=0;i<M;i++){
              U_final[M-2][i] = U_recv[J-2][i];
              U_final[M-1][i] = U_recv[J-1][i];
           }
        }
     }
   ofstream out {"U_horiz_t0.csv"};
   out<<fixed<<setprecision(4);

   for(int i = 0;i<M;i++){
      for(int j=0; j<M;j++){
         out<< U_final[i][j]<<" ";
      }
      out<<endl;
   }
   out.close();
   }

     for (int t=1; t<=N; ++t){
      //Update middle values.
      for (int i=1; i<J-1; ++i){
         for (int j=1; j<M-1; ++j){		
            U[2][i][j] = 2*U[1][i][j] - U[0][i][j]	+ dtdy2*( U[1][i+1][j] + U[1][i-1][j] + U[1][i][j+1] + U[1][i][j-1] - 4*U[1][i][j] ); 	
            }		
         }

      // Swaps top to bottom.
      for(int i=0; i<M-2;i++){
         to_send[i] = U[2][J-2][i+1];
      }
      if(rank != nproc-1){
         MPI_Bsend(to_send,M-2,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD);
         //cout<<"rank "<<rank<<" sent top to bottom fine."<<endl;
      }
      if(rank != 0){
         MPI_Recv(to_recv,M-2,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
         //cout<<"rank "<<rank<<" received top to bottom fine."<<endl;
         for(int i = 0;i<M-2;i++){
            U[2][0][i+1] = to_recv[i];
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);

      // Swaps bottom to top.
      for(int i=0; i<M-2;i++){
         to_send[i] = U[2][1][i+1];
      }
      if(rank != 0){
         MPI_Bsend(to_send,M-2,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD);
         //cout<<"rank "<<rank<<" sent bottom to top fine."<<endl;
      }
      if(rank != nproc - 1){
         MPI_Recv(to_recv,M-2,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
         //cout<<"rank "<<rank<<" received bottom to top fine."<<endl;
         for(int i = 0;i<M-2;i++){
            U[2][J-1][i+1] = to_recv[i];
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);           
      // update the previous times.
      for (int i=0; i<J; ++i){
         for (int j=0; j<M; ++j){		
               U[0][i][j] = U[1][i][j];
               U[1][i][j] = U[2][i][j];
            }		
         }

      if(t==t1 || t==t2 || t==t3){	
         stringstream ss;
         ss << fixed << setprecision(2) << t*dt; // this ensures that the double value gets converted
         string time = ss.str();						 // to string with only 2 trailing digits.

         // Printing initial conditions 
      if(rank !=0){
         MPI_Ssend(U,M*J,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
      }

      if(rank == 0){
         for(int i = 0;i<J-2;i++){
            for(int j =0; j<M;j++){
               U_final[i][j] = U[0][i][j];
            }
         }
         for(int send_rank = 1;send_rank<nproc;send_rank++){
            MPI_Recv(U_recv,M*J,MPI_DOUBLE,send_rank,send_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            
            for(int i =0; i<J-2;i++){
               for(int j=0;j<M;j++){
                  U_final[send_rank*(J-2)+i][j] = U_recv[i][j];
               }
            }
            if(send_rank == nproc-1){
               for(int i=0;i<M;i++){
                  U_final[M-2][i] = U_recv[J-2][i];
                  U_final[M-1][i] = U_recv[J-1][i];
               }
            }
         }
         ofstream out {"U_horiz_t"+ss.str()+".csv"};
         out<<fixed<<setprecision(4);

         for(int i = 0;i<M;i++){
            for(int j=0; j<M;j++){
               out<< U_final[i][j]<<" ";
            }
            out<<endl;
         }
         out.close();
         }
      }
      //cout<<"Iteration "<<t<<" done."<<endl;
   }   

   MPI_Barrier(MPI_COMM_WORLD);
   double t_end = MPI_Wtime();
   double run_time = t_end-t_start;
   cout<<run_time<<endl;
   MPI_Finalize();
}
