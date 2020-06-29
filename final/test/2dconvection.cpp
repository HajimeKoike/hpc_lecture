#include <iostream>
#include <cmath>

using namespace std;

int main(){
	const int nx = 41;
	const int ny = 41;
	const int nt = 25;
	const double c = 1.;
	const double dx = 2./(nx-1);
	const double dy = 2./(ny-1);
	const double dt = dx * .2;

	int i,j,t;

	double u[ny][nx],u_old[ny][nx];

	for ( j=0; j<ny; j++ )
		for ( i=0; i<nx; i++ )
			u[j][i] = 1.;

	for ( j=int(.5/dy); j<int(1/dy+1); j++ )
		for ( i=int(.5/dx); i<int(1/dx+1); i++ )
			u[j][i] = 2.;

	for ( t=0; t<nt; t++ )
	{
		for ( j=0; j<ny; j++ )
		{
			for ( i=0; i<nx; i++ )
			{
				u_old[j][i] = u[j][i];
				cout<<u_old[j][i]<<",";
			}
			cout<<endl;
		}
		for ( j=0; j<ny; j++ )
			for ( i=0; i<nx; i++ )
				u[j][i] = u_old[j][i] - c * dt/dx * ( u_old[j][i] - u_old[j][(i+nx-1)%nx] )
									- c * dt/dy * ( u_old[j][i] - u_old[(j+ny-1)%ny][i] );
	}
}
