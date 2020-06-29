#include <iostream>
#include <cmath>

using namespace std;

int main(){
	const int nx = 11;
	const int ny = 11;
	const int nt = 100;
	const double nu = .01;
	const double dx = 2./(nx-1);
	const double dy = 2./(ny-1);
	const double dt = dx * dy / nu * .25;

	int i,j,t;

	double p[ny][nx],p_old[ny][nx];

	for ( j=0; j<ny; j++ )
	{
		for ( i=0; i<nx; i++ )
		{
			p[j][i] = 0.;
		}
	}

	for ( j=0; j<ny; j++ )
	{
		p[j][nx-1] = j;
	}

	for ( i=0; i<nx; i++ )
	{
		p[0][i] = p[1][i];
		p[ny-1][i] = p[ny-2][i];
	}



	for ( t=0; t<nt; t++ )
	{
		for ( j=0; j<ny; j++ )
		{
			for ( i=0; i<nx; i++ )
			{
				p_old[j][i] = p[j][i];
				cout<<p_old[j][i]<<",";
			}
			cout<<endl;
		}
		for ( j=0; j<ny; j++ )
		{
			for ( i=0; i<nx; i++ )
			{
				p[j][i] =   ( (dy * dy) * ( p_old[j][(i+1)%nx] + p_old[j][(i+nx-1)%nx] )
			 						  + (dx * dx) * ( p_old[(j+1)%ny][i] + p_old[(j+ny-1)%ny][i] ) ) / (2 * (dx*dx + dy*dy));
			}
		}
		for ( j=0; j<ny; j++ )
		{
			p[j][nx-1] = j;
		}

		for ( i=0; i<nx; i++ )
		{
			p[0][i] = p[1][i];
			p[ny-1][i] = p[ny-2][i];
		}

	}
}
