#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 8;
  double x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
	__m256 xvec = _mm256_load_ps(x);
	__m256 yvec = _mm256_load_ps(y);
	__m256 xi = _mm256_set1_ps(x[i]);
	__m256 yi = _mm256_set1_ps(y[i]);
	__m256 mvec = _mm256_load_ps(m);
	__m256 xmask = _mm256_cmp_ps(xvec,xi,_CMP_GT_OQ);
	__m256 ymask = _mm256_cmp_ps(yvec,yi,_CMP_GT_OQ);
	__m256 mask = _mm256_mul_ps(xmask,ymask);	
	xvec = _mm256_blendv_ps(xvec,xvec,mask);
	yvec = _mm256_blendv_ps(yvec,yvec,mask);
	__m256 rxvec = _mm256_sub_ps(xi,xvec);
	__m256 ryvec = _mm256_sub_ps(yi,yvec);
	__m256 rvec  = _mm256_add_ps(_mm256_mul_ps(rxvec,rxvec),_mm256_mul_ps(ryvec,ryvec));
	__m256 r_sqrtrec = _mm256_rsqrt_ps(rvec);
	__m256 fxvec = _mm256_load_ps(fx);
	__m256 fyvec = _mm256_load_ps(fy);
	__m256 r_sqrtrec_cube = _mm256_mul_ps(r_sqrtrec,_mm256_mul_ps(r_sqrtrec,r_sqrtrec));
	fxvec = _mm256_sub_ps(fxvec,_mm256_mul_ps(mvec,r_sqrtrec_cube));
	fyvec = _mm256_sub_ps(fyvec,_mm256_mul_ps(mvec,r_sqrtrec_cube));
	_mm256_store_ps(fx,fxvec);
	_mm256_store_ps(fy,fyvec);
    }
  for(int i=0;i<N;i++){
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
