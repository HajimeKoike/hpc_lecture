#include <iostream>
#include <cmath>

using namespace std;

int main(){
	const int nx = 41;
	const int ny = 41;
	const int nt = 25;
	const double nu = .01;
	const double dx = 2./(nx-1);
	const double dy = 2./(ny-1);
	const double dt = dx * .2;

	int i,j,t;

	double u[ny][nx],u_old[ny][nx];
	double v[ny][nx],v_old[ny][nx];

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
				u[j][i] = u_old[j][i] - u_old[j][i] * dt/dx * ( u_old[j][i] - u_old[j][(i+nx-1)%nx] )
									- v_old[j][i] * dt/dy * ( u_old[j][i] - u_old[(j+ny-1)%ny][i] )
									+ nu * dt/(dx*dx) * (u_old[j][(i+1)%nx] - 2 * u_old[j][i] + u_old[j][(i+nx-1)%nx] )
									+ nu * dt/(dy*dy) * (u_old[(j+1)%ny][i] - 2 * u_old[j][i] + u_old[(j+ny-1)%ny][i] );
				v[j][i] = v_old[j][i] - u_old[j][i] * dt/dx * ( v_old[j][i] - v_old[j][(i+nx-1)%nx] )
									- v_old[j][i] * dt/dy * ( v_old[j][i] - v_old[(j+ny-1)%ny][i] )
									+ nu * dt/(dx*dx) * (v_old[j][(i+1)%nx] - 2 * v_old[j][i] + v_old[j][(i+nx-1)%nx] )
									+ nu * dt/(dy*dy) * (v_old[(j+1)%ny][i] - 2 * v_old[j][i] + v_old[(j+ny-1)%ny][i] );
									
	}
}
