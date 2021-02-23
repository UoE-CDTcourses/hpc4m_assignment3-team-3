#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]){

int rank, ierr, nproc;
MPI_Comm comm;
const int root = 0;
comm = MPI_COMM_WORLD;
MPI_Init(NULL,NULL);
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &nproc);

MPI_Barrier(comm);
// Start timing the code
double t_start = MPI_Wtime();

int M = 101;//2306;  // M length intervals.
int J = (M-1)/nproc + 2; // Length of interval for each process
double T = 1;  // final time.
double dt = 0.2/M;
int N = T/dt;
int t1=0.333/dt, t2=0.666/dt, t3=N; // points at which we want to print our results.
double dy = 2./M;
double dy2 = dy*dy;
double dtdy2 = dt*dt/dy2;

double U_sol[M+1][M+1];
double U[3][M+1][J];
double U_tosend[M+1];
double U_torec[M+1];
double U_tosendup[M+1];
double U_torecup[M+1];

// initialize numerical array with given conditions.
for (int i=0; i<=M; ++i){
	for (int j=0; j<=J; ++j){
		U[0][i][j] = U[1][i][j] = exp( -40 * ( (i*dy-1-0.4)*(i*dy-1-0.4) + ((j+rank*(J-2))*dy-1)*((j+rank*(J-2))*dy-1) ) );
	}
}

for (int i=0; i<J; ++i){ 
	U[0][0][i] = U[0][M][0] = U[1][0][i] = U[1][M][i] = U[2][0][i] = U[2][M][i] = 0; 	
}
if (rank == root)
{
	for (int i=0; i<=M; ++i){
	U[0][i][0]  = U[1][i][0] = U[2][i][0] = 0;
	}
}
if (rank == nproc - 1)
{
	for (int i=0; i<=M; ++i){
	U[0][i][J-1]  = U[1][i][J-1] = U[2][i][J-1] = 0;
	}
}
// changed up to here
//Sending intial arrays back to root process
if (rank != root)
{
	//cout << "Send rank " << rank << endl;
	MPI_Send(&U[0][0][0],J*(M+1),MPI_DOUBLE,root,0,MPI_COMM_WORLD);
	//cout << "Send rank " << rank << endl;
}
// Root process collects all initial arrays and prints them
if (rank == root)
{
	for(int i = 0; i < M+1; i++)
	{
		for(int j = 0; j <J; j++)
		{
			U_sol[i][j] = U[0][i][j];
		}
	}
	for (int i=1; i<nproc; i++)
	{
		int k = 0;
		MPI_Recv(&U_torec,J*(M+1),MPI_DOUBLE,i,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(int row = 0; row < M+1; row++)
		{
			for(int col = 0; col < J ;col++)
			{
				U_sol[row][(J-2) * i + col] = U_torec[k];
				cout << U_torec[k] << endl;
				k++;
			}
		}
	}

	//print initial U values to file, row by row.
	ofstream out {"U_t0ser.csv"};
	out<<fixed<<setprecision(4);
		for(int i=0; i<=M; ++i){
			for(int j=0; j<=M; ++j){					
				out<<U_sol[i][j]<<" ";
				//cout << "Here?\n i = "<< i << " j = " << j <<"\t" << U_sol[i][j]<< endl;
			}
			out<<endl;
		}	//cout << "HEre? rank -" << rank << endl;
	out.close();
}


// use numerical scheme to obtain the future values of U.
for (int t=1; t<=N; ++t){
	
	for (int i=1; i<M; ++i){
		for (int j=1; j<J; ++j){		
			U[2][i][j] = 2*U[1][i][j] - U[0][i][j]	 
						 	 + dtdy2*( U[1][i+1][j] + U[1][i-1][j] + U[1][i][j+1] + U[1][i][j-1] - 4*U[1][i][j] ); 	// this has been changed
		}		
	}
	if (rank != root){    // every procsses except root sends all elements in boundary to previous process...
			//MPI_Send(&U[2][1][0],M+1, MPI_DOUBLE, rank-1, 0, comm);
			for(int i = 0; i < M+1;i++)
			{
				U_tosend[i] = U[2][i][1];
			}
			MPI_Send(&U_tosend,M+1, MPI_DOUBLE, rank - 1, 0, comm);
	}

		if (rank != nproc-1){  // ...every processes except the final processes receives these points.
			//MPI_Recv(&U[2][J-1][0],M+1, MPI_DOUBLE, rank+1, 0, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&U_torec,M+1, MPI_DOUBLE, rank+1, 0, comm, MPI_STATUS_IGNORE);
			for(int i = 0; i < M+1;i++)
			{
				U[2][i][J-1] = U_torec[i];
			}
		}
	if (rank !=nproc-1){    // every procsses except the final process sends all elements in boundary to next process...
			//MPI_Send(&U[2][J-2][0],M+1, MPI_DOUBLE, rank+1, 0, comm);
			for(int i = 0; i < M+1;i++)
			{
				U_tosendup[i] = U[2][i][J-2];
			}
			MPI_Send(&U_tosendup,M+1, MPI_DOUBLE, rank+1, 0, comm);
	}

		if (rank != root){  // ...every processes except the root processes receives these points.
			//MPI_Recv(&U[2][0][0],M+1, MPI_DOUBLE, rank-1, 0, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&U_torecup,M+1, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE);
			for(int i = 0; i < M+1;i++)
			{
				U[2][i][0] = U_torecup[i];
			}
		}
	
	// update the previous times.
	for (int i=0; i<M+1; ++i){
		for (int j=0; j<J; ++j){		
			U[0][i][j] = U[1][i][j];
			U[1][i][j] = U[2][i][j];
			}		
	}

if(t==t1 || t==t2 || t==t3){	
if (rank != root)
{
	MPI_Send(&U[2][0][0],J*(M+1),MPI_DOUBLE,root,0,MPI_COMM_WORLD);
}
// Root process collects all initial arrays and prints them
if (rank == root)
{
	for(int i = 0; i < M+1; i++)
	{
		for(int j = 0; j <J; j++)
		{
			U_sol[i][j] = U[2][i][j];
		}
	}
	for (int i=1; i<nproc; i++)
	{
		int k = 0;
		MPI_Recv(&U_torec,J*(M+1),MPI_DOUBLE,i,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(int row = 0; row < M+1; row++)
		{
			for(int col = 0; col < J ;col++)
			{
				U_sol[row][(J-2) * i + col] = U_torec[k];
				k++;
			}
		}
	}	
   // print out files for fixed times	
		stringstream ss;
		ss << fixed << setprecision(2) << t*dt; // this ensures that the double value gets converted
		string time = ss.str();						 // to string with only 2 trailing digits.
		
		ofstream out {"U_t"+ss.str()+"ser.csv"};
		out<<fixed<<setprecision(4);			
			for(int i=0; i<=M; ++i){
				for(int j=0; j<=M; ++j){		
					out<<U_sol[i][j]<<" ";
				}
				out<<endl;
			}
		out.close();
	}
	cout<<"iteration "<< t <<" done"<<endl;	

}}

//Ensure all processes have reached this point
MPI_Barrier(comm);
// Record time taken for the code to run
double t_end = MPI_Wtime();
double run_time = t_end - t_start;

if (rank == root)
{
	cout<< "\nThe code took  "<< run_time << " seconds to run.\n";
}

return 0;
MPI_Finalize();
}