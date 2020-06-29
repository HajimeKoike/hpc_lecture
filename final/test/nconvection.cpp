#include <iostream>
#include <cmath>

using namespace std;

int main(){
	const int nx = 41;
	const double dx = 2./(nx-1);
	const int nt = 25;
	const double dt = .025;
	const double c = 1.;

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
			u[i+1] = u_old[i] - u_old[i] * dt/dx * ( u_old[i] - u_old[(i+nx-1)%nx] );

		cout<<endl;
	}
}
