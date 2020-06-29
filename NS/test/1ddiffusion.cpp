#include <iostream>
#include <cmath>

using namespace std;

int main(){
	const int nx = 41;
	const double dx = 2./(nx-1);
	const int nt = 25;
	const double nu = 1.;
	const double sigma = .2;
	const double dt = sigma * dx * dx / nu;

	int i,t;

	double u[nx],u_old[nx];

	for ( i=0; i<nx; i++ )
		u[i] = 1.;
	u[int(.5/dx)] = 2.;

	for ( t=0; t<nt; t++ )
	{
		for ( i=0; i<nx; i++ )
		{
			u_old[i] = u[i];
			cout<<u_old[i]<<",";
		}

		for ( i=0; i<nx; i++ )
			u[(i+1)%nx] = u_old[i] + nu * dt / (dx * dx) * ( u_old[(i+1)%nx] - 2 * u_old[i] + u_old[(i+nx-1)%nx] );

		cout<<endl;
	}
}
