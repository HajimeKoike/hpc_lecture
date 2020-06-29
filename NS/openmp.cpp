#include <iostream>
#include <cmath>

using namespace std;

int main(){
	const int nx = 41;
	const int ny = 41;
	const int nt = 10;
	const int nit = 50;
	const double c = 1.;
	const double rho = 1.;
	const double nu = .1;
	const double F = 1.;
	const double dx = 2./(nx-1);
	const double dy = 2./(ny-1);
	const double dt = .01;

	int i,j,t,q;

	double u[ny][nx],u_old[ny][nx];
	double v[ny][nx],v_old[ny][nx];
	double p[ny][nx],p_old[ny][nx];
	double b[ny][nx];

	double u_sum,u_old_sum,u_diff;

	for ( j=0; j<ny; j++ )
	{
#pragma omp parallel for
		for ( i=0; i<nx; i++ )
		{
			u[j][i] = 0.;
			v[j][i] = 0.;
			p[j][i] = 1.;
			b[j][i] = 0.;
		}
	}

	u_diff = 1.;
	t = 0;
	for ( t=0; t<nt; t++)
	{
		for ( j=0; j<ny; j++)
		{
#pragma omp parallel for
			for ( i=0; i<nx; i++)
			{
				u_old[j][i] = u[j][i];
				v_old[j][i] = v[j][i];
			}
		}

		//Source
		for ( j=0; j<ny; j++ )
		{
#pragma omp parallel for
			for ( i=0; i<nx; i++ )
			{
				b[j][i] = ( rho * (1 /dt * (u[j][(i+1)%nx] - u[j][(i+nx-1)%nx] ) / (2*dx)
							+ ( v[(j+1)%ny][i] - v[(j+ny-1)%ny][i] ) / (2*dy) )
						- pow( ( u[j][(i+1)%nx] - u[j][(i+nx-1)%nx] ) / (2*dx), 2)
						- 2 * ( ( u[(j+1)%ny][i] - u[(i+ny-1)%ny][i] ) / (2*dy) 
									* ( v[j][(i+1)%nx] - v[j][(i+nx-1)%nx] ) / (2*dx) ) 
						- pow( ( v[(j+1)%ny][i]-v[(j+ny-1)%ny][i] ) / (2*dy), 2) );
				cout<<b[j][i]<<",";
			}
			cout<<endl;
		}


		//Pressure
		for ( q=0; q<nit; q++ )
		{
			for ( j=0; j<ny; j++ )
			{
#pragma omp parallel for
				for ( i=0; i<nx; i++ )
				{
					p_old[j][i] = p[j][i];
					cout<<p_old[j][i]<<",";
				}
				cout<<endl;
			}
			for ( j=0; j<ny; j++ )
			{
#pragma omp parallel for
				for ( i=0; i<nx; i++ )
				{
					p[j][i] =   ( (dy * dy) * ( p_old[j][(i+1)%nx] + p_old[j][(i+nx-1)%nx] )
										  + (dx * dx) * ( p_old[(j+1)%ny][i] + p_old[(j+ny-1)%ny][i] ) ) / (2 * (dx*dx + dy*dy)) 
								- b[j][i] * dx*dx * dy*dy  / (2 * (dx*dx + dy*dy));
				}
			}
			//
#pragma omp parallel for
 			for ( j=0; j<ny; j++ )
			{
				p[j][nx-1] = p[j][nx-2];
				p[j][0] = p[j][1];
			}

#pragma omp parallel for
			for ( i=0; i<nx; i++ )
			{
				p[0][i] = p[1][i];
				p[ny-1][i] = p[ny-2][j];
			}
		}

		//Cavity
		for ( j=0; j<ny; j++)
		{
#pragma omp parallel for
			for ( i=0; i<nx; i++ )
			{
			u[j][i] = ( u_old[j][i] - u_old[j][i] * dt/dx * ( u_old[j][i] - u_old[j][(i+nx-1)%nx] ) 
									- v_old[j][i] * dt/dy * ( u_old[j][i] - u_old[(j+ny-1)%ny][i] )
					- dt / (2 * rho * dx) * ( p[j][(i+1)%nx] - p[j][(j+nx-1)%nx] ) 
					+ nu * (dt / dx*dx * ( u_old[j][(i+1)%nx] - 2 * u_old[j][i] + u_old[j][(i+nx-1)%nx] ) 
						+ dt / (dy*dy) * ( u_old[(j+1)%ny][i] - 2 * u_old[j][i] + u_old[(j+ny-1)%ny][i]) ) 
					+ F * dt );

			v[j][i] = ( v_old[j][i] - u_old[j][i] * dt/dx * ( v_old[j][i] - v_old[j][(i+nx-1)%nx] ) 
									- v_old[j][i] * dt/dy * ( v_old[j][i] - v_old[(j+ny-1)%ny][i] )
					- dt / (2 * rho * dy) * ( p[(j+1)%ny][i] - p[(j+ny-1)%ny][j] ) 
					+ nu * (dt / dx*dx * ( v_old[j][(i+1)%nx] - 2 * v_old[j][i] + v_old[j][(i+nx-1)%nx] ) 
						+ dt / (dy*dy) * ( v_old[(j+1)%ny][i] - 2 * v_old[j][i] + v_old[(j+ny-1)%ny][i]) ) );
			}
		}

#pragma omp parallel for
		for ( i=0; i<nx; i++ )
		{
			u[0][i] = 0.;
			u[ny-1][i] = 0.;
			v[0][i] = 0.;
			v[ny-1][i] = 0.;
		}

#pragma omp parallel for
		for ( j=0; j<ny; j++)
		{
			u[j][0] = 0.;
			u[j][nx-1] = 0.;
			v[j][0] = 0.;
			v[j][nx-1] = 0.;
		}

		u_sum = 0.;
		u_old_sum = 0.;
		u_diff = 0.;
		for ( j=0; j<ny; j++ )
		{
#pragma omp parallel for reduction(+:u_sum)
#pragma omp parallel for reduction(+:u_old_sum)
			for ( i=0; i<nx; i++ )
			{
				u_sum += u[j][i];
				u_old_sum += u_old[j][i];
			}
		}
		u_diff = (u_sum - u_old_sum) / u_sum;
		t++;

	}

}
