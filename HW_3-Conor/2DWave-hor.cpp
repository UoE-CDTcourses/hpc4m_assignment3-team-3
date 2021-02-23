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

int M = 2301;//2306;  // M length intervals.
int J = (M-1)/nproc + 2; // Length of interval for each process
double T = 1;  // final time.
double dt = 0.2/M;
int N = T/dt;
int t1=0.333/dt, t2=0.666/dt, t3=N; // points at which we want to print our results.
double dy = 2./M;
double dy2 = dy*dy;
double dtdy2 = dt*dt/dy2;

double** U_sol = new double* [M+1]; // this creates an array of size 3 which is ready do hold 2D arrays in each entry.
for (int i = 0; i < M+1; ++i) {   // this loop puts an array of size M+1 on each of the three entries.
  U_sol[i] = new double[M+1];
  }


// For vertical strips make size U[3][J][M+1]
double*** U = new double** [3]; // this creates an array of size 3 which is ready do hold 2D arrays in each entry.
for (int i = 0; i < 3; ++i) {   // this loop puts an array of size M+1 on each of the three entries.
  U[i] = new double*[J];		  // we get a 2D array.
  for (int j = 0; j < M+1; ++j){ // put another array on each entry of the 2D array to get a 3D array.
    U[i][j] = new double[M+1];
  }
}

// initialize numerical array with given conditions.
for (int i=0; i<J; ++i){
	for (int j=0; j<=M; ++j){
		U[0][i][j] = U[1][i][j] = exp( -40 * ( ((i+rank*(J-2))*dy-1-0.4)*((i+rank*(J-2))*dy-1-0.4) + (j*dy-1)*(j*dy-1) ) );
	}
}

for (int i=0; i<J; ++i){ 
	U[0][i][0] = U[0][i][M] = U[1][i][0] = U[1][i][M] = U[2][i][0] = U[2][i][M] = 0; 	
}
if (rank == root)
{
	for (int i=0; i<=M; ++i){
	U[0][0][i]  = U[1][0][i] = U[2][0][i] = 0;
	}
}
if (rank == nproc - 1)
{
	for (int i=0; i<=M; ++i){
	U[0][J-1][i]  = U[1][J-1][i] = U[2][J-1][i] = 0;
	}
}

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
	for(int i = 0; i < J; i++)
	{
		for(int j = 0; j <=M; j++)
		{
			U_sol[i][j] = U[0][i][j];
			//cout << "I haven't made it yet - Rank = " << rank << " i = " << i << "j = " << j<< endl;
		}
	}
	for (int i=1; i<nproc; i++)
	{
		MPI_Recv(&U_sol[(J-2) * i][0],J*(M+1),MPI_DOUBLE,i,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//print initial U values to file, row by row.
	ofstream out {"U_t0ser.csv"};
	out<<fixed<<setprecision(4);
		for(int i=0; i<=M; ++i){
			for(int j=0; j<=M; ++j){					
				out<<U_sol[i][j]<<" ";
			}
			out<<endl;
		}	
	out.close();
}


// use numerical scheme to obtain the future values of U.
for (int t=1; t<=N; ++t){
	
	for (int i=1; i<J; ++i){
		for (int j=1; j<M; ++j){		
			U[2][i][j] = 2*U[1][i][j] - U[0][i][j]	 
						 	 + dtdy2*( U[1][i+1][j] + U[1][i-1][j] + U[1][i][j+1] + U[1][i][j-1] - 4*U[1][i][j] ); 	
		}		
	}
	if (rank != root){    // every procsses except root sends all elements in boundary to previous process...
			MPI_Send(&U[2][1][0],M+1, MPI_DOUBLE, rank-1, 0, comm);
	}

		if (rank != nproc-1){  // ...every processes except the final processes receives these points.
			MPI_Recv(&U[2][J-1][0],M+1, MPI_DOUBLE, rank+1, 0, comm, MPI_STATUS_IGNORE);
		}

	if (rank !=nproc-1){    // every procsses except the final process sends all elements in boundary to next process...
			MPI_Send(&U[2][J-2][0],M+1, MPI_DOUBLE, rank+1, 0, comm);
	}

		if (rank != root){  // ...every processes except the root processes receives these points.
			MPI_Recv(&U[2][0][0],M+1, MPI_DOUBLE, rank-1, 0, comm, MPI_STATUS_IGNORE);
		}
	
	// update the previous times.
	for (int i=0; i<J; ++i){
		for (int j=0; j<M+1; ++j){		
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
	for(int i = 0; i < J; i++)
	{
		for(int j = 0; j <=M; j++)
		{
			U_sol[i][j] = U[2][i][j];
		}
	}
	for (int i=1; i<nproc; i++)
	{
		MPI_Recv(&U_sol[(J-2) * i][0],J*(M+1),MPI_DOUBLE,i,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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


// now we need to delete the U_sol array.
MPI_Barrier(comm);
for (int i = 0; i < M+1; i++)
{
    delete[] U_sol[i];
}
delete[] U_sol;
// now we need to delete the U array.
for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < J; j++){
	 }  
    delete[] U[i];
	MPI_Barrier(comm);
}
delete[] U;
return 0;
MPI_Finalize();
}