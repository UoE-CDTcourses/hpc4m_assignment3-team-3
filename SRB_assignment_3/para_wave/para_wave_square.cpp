#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

static const double pi = 3.1415926536;
static const int root_nproc = 2;
static const int M = 2306;
static const double T = 1;

int main(int argc, char* argv[]){

   int rank, nproc, ierr, rank1, rank2;  
   double dt = 0.2/(M);
   int N = T/dt;
   int t1=0.333/dt;
   int t2=0.666/dt, t3=N; 
   double dy = 2./(M-1);
   double dx = dy;
   double dy2 = dy*dy;
   double dtdy2 = dt*dt/dy2;
   double x,y;
   int column,row;
   MPI_Comm comm;

   MPI_Init(NULL,NULL);

   MPI_Barrier(MPI_COMM_WORLD);
   double t_start = MPI_Wtime();

   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);


   int J = 2+(M-2)/root_nproc;
   double to_send[J-1],to_recv[J-1];

   double U[3][J][J];
   double U_final[M][M];
   double U_recv[J][J];
   int sender_row;
   int sender_column;

   // Which row/column is the square region corresponding to the given process located on?
   row = floor(rank/root_nproc);
   column = rank % root_nproc; 

   // Functional initial condition.
       for (int i=0; i<J; ++i){
         for (int j=0; j<J; ++j){

            y = (row*(J-2)+i)*dy - 1;
            x = (column*(J-2)+j)*dx-1;

            U[0][i][j] = U[1][i][j] = exp( -40 *  ((x-0.4)*(x-0.4) + (y*y )));
         }
      } 
   // Top boundaries required for upper row of squares.
   if(row == 0){
      for(int i=0;i<J;i++){
         U[0][0][i] = 0;
         U[1][0][i] = 0;
         U[2][0][i] = 0;
      } 
   }

   // Bottom boundaries required for lower row of squares.
   if(row == root_nproc-1){
      for(int i=0;i<J;i++){
         U[0][J-1][i] = 0;
         U[1][J-1][i] = 0;
         U[2][J-1][i] = 0;
      } 
   }

   // Left boundaries required for left column of squares.
   if(column == 0){
      for(int i=0;i<J;i++){
         U[0][i][0] = 0;
         U[1][i][0] = 0;
         U[2][i][0] = 0;
      } 
   }

   // Right boundaries required for right column of squares.
   if(column == sqrt(nproc)-1){
      for(int i=0;i<J;i++){
         U[0][i][J-1] = 0;
         U[1][i][J-1] = 0;
         U[2][i][J-1] = 0;
      } 
   }
   // Print out initial conditions to file.
   if(rank != 0){  
      cout<<"rank"<<rank<<"sending initial condition"<<endl;
      MPI_Ssend(U,J*J,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
   }  

   if(rank == 0){
      for(int i=0;i<J-2;i++){
         for(int j=0;j<J-2;j++){
            U_final[i][j] = U[0][i][j];
         }
      }
      for(int send_rank = 1;send_rank<nproc;send_rank++){
         MPI_Recv(U_recv,J*J,MPI_DOUBLE,send_rank,send_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
         cout<<"Received rank "<<send_rank<<" initial condition"<<endl; 

         sender_row = floor(send_rank/root_nproc);
         sender_column = send_rank % root_nproc;

         for(int i =0; i<J-2;i++){
            for(int j=0; j<J-2;j++){
               U_final[sender_row*(J-2)+i][sender_column*(J-2)+j] = U_recv[i][j];
            }
         }

         if(sender_row == root_nproc-1){
            for(int i = 0; i<J-2;i++){
               U_final[M-2][sender_column*(J-2)+i] = U_recv[J-2][i];
               U_final[M-1][sender_column*(J-2)+i] = U_recv[J-1][i];
            }
         }

         if(sender_column == root_nproc-1){
            for(int i = 0;i<J-2;i++){
               U_final[sender_row*(J-2)+i][M-2] = U_recv[i][J-2];
               U_final[sender_row*(J-2)+i][M-1] = U_recv[i][J-1];
            }
         }

         if(send_rank == nproc - 1){
            U_final[M-1][M-2] = U_recv[J-1][J-2];
            U_final[M-2][M-1] = U_recv[J-2][J-1];
            U_final[M-1][M-1] = U_recv[J-1][J-1];
         }
      }

      ofstream out {"U_square_t0.csv"};
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
         for (int j=1; j<J-1; ++j){		
            U[2][i][j] = 2*U[1][i][j] - U[0][i][j]	+ dtdy2*( U[1][i+1][j] + U[1][i-1][j] + U[1][i][j+1] + U[1][i][j-1] - 4*U[1][i][j] ); 	
            }		
         }
      MPI_Barrier(MPI_COMM_WORLD);
      // Halo swapping.
         // Swaps left to right.
         for(int i=0; i<J-2;i++){
            to_send[i] = U[2][i+1][J-2];
         }
         if(column != root_nproc-1){
            MPI_Bsend(to_send,J-2,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD);
            //cout<<"rank "<<rank<<" sent left to right fine."<<endl;
         }
         if(column != 0){
            MPI_Recv(to_recv,J-2,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout<<"rank "<<rank<<" received left to right fine."<<endl;
            for(int i = 0;i<J-2;i++){
               U[2][i+1][0] = to_recv[i];
            }
         }
         MPI_Barrier(MPI_COMM_WORLD);

         // Swaps right to left.
         for(int i=0; i<J-2;i++){
            to_send[i] = U[2][i+1][1];
         }
         if(column != 0){
            MPI_Bsend(to_send,J-2,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD);
            //cout<<"rank "<<rank<<" sent right to left fine."<<endl;
         }
         if(column != root_nproc-1){
            MPI_Recv(to_recv,J-2,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout<<"rank "<<rank<<" received right to left fine."<<endl;
            for(int i = 0;i<J-2;i++){
               U[2][i+1][J-1] = to_recv[i];
            }
         }
         MPI_Barrier(MPI_COMM_WORLD);

         // Swaps top to bottom.
         for(int i=0; i<J-1;i++){
            to_send[i] = U[2][J-2][i+1];
         }
         if(row != root_nproc-1){
            MPI_Bsend(to_send,J-1,MPI_DOUBLE,rank+root_nproc,rank,MPI_COMM_WORLD);
            //cout<<"rank "<<rank<<" sent top to bottom fine."<<endl;
         }
         if(row != 0){
            MPI_Recv(to_recv,J-1,MPI_DOUBLE,rank-root_nproc,rank-root_nproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout<<"rank "<<rank<<" received top to bottom fine."<<endl;
            for(int i = 0;i<J-1;i++){
               U[2][0][i+1] = to_recv[i];
            }
         }
         MPI_Barrier(MPI_COMM_WORLD);

         // Swaps bottom to top.
         for(int i=0; i<J-1;i++){
            to_send[i] = U[2][1][i+1];
         }
         if(row != 0){
            MPI_Bsend(to_send,J-1,MPI_DOUBLE,rank-root_nproc,rank,MPI_COMM_WORLD);
            //cout<<"rank "<<rank<<" sent bottom to top fine."<<endl;
         }
         if(row != root_nproc - 1){
            MPI_Recv(to_recv,J-1,MPI_DOUBLE,rank+root_nproc,rank+root_nproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout<<"rank "<<rank<<" received bottom to top fine."<<endl;
            for(int i = 0;i<J-1;i++){
               U[2][J-1][i+1] = to_recv[i];
            }
         }
         MPI_Barrier(MPI_COMM_WORLD);

      // update the previous times.
      for (int i=0; i<J; ++i){
         for (int j=0; j<J; ++j){		
               U[0][i][j] = U[1][i][j];
               U[1][i][j] = U[2][i][j];
            }		
         }
      //cout<<"rank "<<rank<<" updated previous timesteps."<<endl;  

      if(t==t1 || t==t2 || t==t3){	
         stringstream ss;
         ss << fixed << setprecision(2) << t*dt; // this ensures that the double value gets converted
         string time = ss.str();						 // to string with only 2 trailing digits.

         if(rank != 0){  
            cout<<"rank"<<rank<<"sending t="+ss.str()<<endl;
            MPI_Ssend(U,J*J,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
         }  

         if(rank == 0){
            for(int i=0;i<J-2;i++){
               for(int j=0;j<J-2;j++){
                  U_final[i][j] = U[0][i][j];
               }
            }
         for(int send_rank = 1;send_rank<nproc;send_rank++){
            MPI_Recv(U_recv,J*J,MPI_DOUBLE,send_rank,send_rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            cout<<"Received rank "<<send_rank<<" t="+ss.str()<<endl; 

            sender_row = floor(send_rank/root_nproc);
            sender_column = send_rank % root_nproc;

            for(int i =0; i<J-2;i++){
               for(int j=0; j<J-2;j++){
                  U_final[sender_row*(J-2)+i][sender_column*(J-2)+j] = U_recv[i][j];
               }
            }  

         if(sender_row == root_nproc-1){
            for(int i = 0; i<J-2;i++){
               U_final[M-2][sender_column*(J-2)+i] = U_recv[J-2][i];
               U_final[M-1][sender_column*(J-2)+i] = U_recv[J-1][i];
            }
         }

         if(sender_column == root_nproc-1){
            for(int i = 0;i<J-2;i++){
               U_final[sender_row*(J-2)+i][M-2] = U_recv[i][J-2];
               U_final[sender_row*(J-2)+i][M-1] = U_recv[i][J-1];
            }
         }

         if(send_rank == nproc - 1){
            U_final[M-1][M-2] = U_recv[J-1][J-2];
            U_final[M-2][M-1] = U_recv[J-2][J-1];
            U_final[M-1][M-1] = U_recv[J-1][J-1];
         }
      }

      ofstream out {"U_square_t"+ss.str()+".csv"};
      out<<fixed<<setprecision(4);

      for(int i = 0;i<M;i++){
         for(int j=0; j<M;j++){
               out<< U_final[i][j]<<" ";
            }
            out<<endl;
         }
      out.close();
      }   
      cout<<"Iteration "<<t<<" done."<<endl; 
   }
}      
double t_end = MPI_Wtime();
double run_time = t_end-t_start;
cout<<run_time<<endl;

MPI_Finalize();
return 0;
}
