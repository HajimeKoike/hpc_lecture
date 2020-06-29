#include <iostream>
#include <cmath>

#define naive 0
#define interesting 1

using namespace std;

int main(){

	const int nx = 101;
	const int nt = 100;
	const double dx = 2 * M_PI / ( nx - 1 );
	const double nu = .07;
	const double dt = dx * nu;

	double u[nx], u_old[nx];

	int i,t;
	
#if naive
	//Naive initial condition
	for ( i=0; i<nx; i++ )
		u[i] = 1.;
	u[int(.5/dx)] = 2.;
#endif
#if interesting
	//Interesting initial condition
	double phi[nx],d_phi[nx];
	for ( i=0; i<nx; i++ )
	{
		phi[i] = exp( - i * i / (4 * nu) ) + exp( -(i - 2 * M_PI) * (i - 2 * M_PI) / (4 * nu) );
		d_phi[i] = - 2 * i / (4 * nu) * exp( -i * i / (4 * nu) ) - 2 * ( i - 2 * M_PI) / (4 * nu) * exp( - (i - 2 * M_PI) * (i - 2 * M_PI) / (4 * nu) );
		u[i] = - 2 * nu * d_phi[i] / phi[i]  + 4;
		cout<<d_phi[i]<<",";
	}
#endif
	cout<<endl;
	return 0;


	for ( t=0; t<nt; t++ )
	{
		for ( i=0; i<nx; i++ )
		{
			u_old[i] = u[i];
			cout<<u_old[i]<<",";
		}
		cout<<endl;	

		for ( i=0; i<nx; i++ )
			u[i] = u_old[i] - u_old[i] * dt / dx * ( u_old[i] - u_old[(i+nx-1)%nx] ) 
				+ nu * dt / (dx*dx) * ( u_old[(i+1)%nx] - 2*u_old[i] + u_old[(i+nx-1)%nx] );

		//u[0] = u_old[0] - u_old[0] * dt / dx * ( u_old[0] - u_old[0-1] ) 
		//	+ nu * dt / (dx*dx) * ( u_old[0+1] - 2*u_old[0] + u_old[0-1] );
		//u[-1] = u[0];
	}
}
