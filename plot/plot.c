#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/file.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <math.h>
#include <acb.h>
#include "smalljac.h"
extern int i_poly_parse(long *,int,char *);
extern acb_t *fft_init(long,long);
extern void fft(acb_t *);
extern void ifft(acb_t *);

static int *a,*ftable;
static smalljac_curve_t curve;
static long maxn;

/* invert the polynomial p to precision k */
static void inverse(int *ans,int *p,int k) {
	int i,j,c[32];

	c[0] = 1;
	for (j=1;j<=k;j++) c[j] = 0;
	for (i=0;i<=k;i++) {
		ans[i] = c[0];
		for (j=1;j<=k-i;j++)
			c[j-1] = c[j]-p[j]*ans[i];
	}
}

// floor(log(x)/log(p))
static inline int logp(unsigned long x,unsigned long p) {
	int k;
	for (k=0;x>=p;k++)
		x /= p;
	return k;
}

static struct {
	long p,f[5];
} bad_lfactors[8];

static int recordlpoly(smalljac_Qcurve_t C,unsigned long p,
int good,long aa[],int n,void *arg) {
	int f[64],c[64],i,j,k;
	unsigned long q;

	f[0] = 1; for (j=1;j<64;j++) f[j] = 0;
	if (good) {
		f[1] = aa[0];
		f[2] = aa[1];
		f[3] = f[1]*p;
		f[4] = p*p;
	} else {
		for (i=0;bad_lfactors[i].p;i++)
			if (bad_lfactors[i].p == p) {
				for (j=0;j<=4;j++)
					f[j] = bad_lfactors[i].f[j];
				break;
			}
		if (!bad_lfactors[i].p) {
			fprintf(stderr,"unexpected bad reduction\n");
			exit(1);
		}
	}

	k = logp(maxn,p);
	inverse(c,f,k);
	for (j=1,q=1;j<=k;j++)
		q *= p, a[q] = c[j];
	return 1;
}

static void allocate(long n) {
	int k;
	static long nallocated;

	if (n <= nallocated) return;

	// allocate space for coefficients and factor table
	if (nallocated) a++, ftable++;
	a = (int *)realloc(a,n*sizeof(*a))-1;
	ftable = (int *)realloc(ftable,n*sizeof(*ftable))-1;
	nallocated = n;

  // sieve to compute first factor table */
	for (n=1;n<=nallocated;n++)
		ftable[n] = n;
	for (k=(int)sqrtl((long double)nallocated);k>=2;k--)
		for (n=k+k;n<=nallocated;n+=k)
			ftable[n] = k;
}

#define prec 300
#define cutoff 270 // truncate S at cutoff*x
void compute_coeffs(double xmax) {
	long p,n,q;

	maxn = (long)(cutoff*xmax);
	allocate(maxn);
	a[1] = 1;
	smalljac_Lpolys(curve,1,maxn,0,recordlpoly,(void *)0);
	for (n=2;n<=maxn;n++) {
		p = ftable[n];
		for (q=p;n/q%p==0;q*=p);
		a[n] = a[q] * a[n/q];
	}
}

extern void K0(arb_t,const arb_t);
static void S(arb_t s,const arb_t x) {
	static int init;
	static arb_t c,t,l;
	long n,m;

	if (!init) {
		arb_init(c); arb_init(t);
		arb_init(l);
		init = 1;
	}

	arb_mul_ui(t,x,cutoff,prec);
	m = arf_get_si(arb_midref(t),ARF_RND_FLOOR);

	arb_rsqrt(c,x,prec);
	arb_const_pi(t,prec);
	arb_mul_2exp_si(t,t,2);
	arb_mul(c,c,t,prec);
	arb_zero(l);
	for (n=1;n<=m;n++) {
		arb_sqrt_ui(t,n,prec);
		arb_mul(t,t,c,prec);
		K0(t,t);
		arb_addmul_si(l,t,a[n],prec);
	}
	arb_div(s,l,x,prec);
}

#define A 5
#define sample_points 1024
#define upsample_factor 16
#define big_sample_points (upsample_factor*sample_points)
#define bigA (A*upsample_factor)
#define fpoints (22*bigA)
static acb_t *data,*bigdata;
static int epsilon,rank;
static long hash;
static arb_t scale;
static arb_t *zero;
static int nzeros,mzeros;
static arb_t f[fpoints];

// abs(Gamma_C(1+I*t))^2
void gamma_factor(arb_t g,const arb_t t) {
	static int init;
	static arb_t temp,pi;
	if (!init) {
		arb_init(temp);
		arb_init(pi);
		arb_const_pi(pi,prec);
		init = 1;
	}
	if (arb_contains_zero(t)) {
		arb_mul(temp,pi,pi,prec);
		arb_inv(g,temp,prec);
	} else {
		arb_mul(temp,t,pi,prec);
		arb_sinh(temp,temp,prec);
		arb_mul(temp,temp,pi,prec);
		arb_div(g,t,temp,prec);
	}
}

static void Lambda(arb_t res,const arb_t t) {
	static int init;
	static arb_t c,At,temp;
	int m,absm,sign;

	if (!init) {
		arb_init(c);
		arb_init(At);
		arb_init(temp);
	}
	arb_mul_ui(At,t,A,prec);
	m = arf_get_si(arb_midref(At),ARF_RND_NEAR);
	arb_sub_si(temp,At,m,prec);
	if (arb_contains_zero(temp)) {
		absm = (m < 0) ? -m : m;
		sign = (m < 0) ? epsilon : 1;
		arb_mul_si(res,acb_realref(data[absm]),sign,prec);
		return;
	}
	arb_sin_pi(c,At,prec);
	arb_const_pi(temp,prec);
	arb_div(c,c,temp,prec);

	arb_zero(res);
	for (m=-sample_points/2;m<=sample_points/2;m++) {
		absm = (m < 0) ? -m : m;
		sign = (m < 0) ? epsilon : 1;
		if (m & 1) sign = -sign;
		arb_sub_si(temp,At,m,prec);
		arb_div(temp,c,temp,prec);
		arb_mul(temp,temp,acb_realref(data[absm]),prec);
		if (sign < 0) arb_neg(temp,temp);
		arb_add(res,res,temp,prec);
	}
}

static void Lambda_adjusted(arb_t res,const arb_t t) {
	static int init;
	static arb_t temp,L;
	int i;

	if (!init) {
		arb_init(temp);
		arb_init(L);
		init = 1;
	}
	Lambda(L,t);
	arb_pow_ui(temp,t,rank,prec);
	if (rank & 2) arb_neg(temp,temp);
	arb_div(L,L,temp,prec);
	for (i=0;i<nzeros;i++) {
		arb_div(temp,t,zero[i],prec);
		arb_mul(temp,temp,temp,prec);
		arb_neg(temp,temp); arb_add_ui(temp,temp,1,prec);
		arb_div(L,L,temp,prec);
	}
	arb_set(res,L);
}

static char *buffer;
static long nbuffer;
static int printstr(const char *s) {
	long n;
	if (!buffer || (n=strlen(buffer)+strlen(s)+1) > nbuffer) {
		buffer = (char *)realloc(buffer,n);
		nbuffer = n;
	}
	strcat(buffer,s);
}

static void print(arb_t x) {
	char *s;
	if (arb_contains_zero(x))
		printstr("0");
	else {
		s = arb_get_str(x,16,ARB_STR_NO_RADIUS);
		printstr(s);
		free(s);
	}
}

static void printZ(void) {
	static int init;
	static arb_t t;
	int m,mmax=30*bigA,absm;
	char buf[16];

	if (!init) {
		arb_init(t);
		init = 1;
	}
	printstr("[");
	for (m=-mmax;m<=mmax;m+=bigA/80) {
		absm = (m<0) ? -m : m;
		arb_set_si(t,m);
		arb_div_ui(t,t,bigA,prec);
		printstr("[");

		sprintf(buf,"%s%d.%04d,",m<0?"-":"",
			absm/bigA,(absm*10000/bigA)%10000);
		printstr(buf);
		gamma_factor(t,t);
		arb_div(t,scale,t,prec);
		arb_mul(t,t,acb_realref(bigdata[absm]),prec);
		if (m < 0 && epsilon < 0)
			arb_neg(t,t);
		if (!m && rank > 0) printstr("0"); else print(t);
		printstr("]");
		printstr(m<mmax ? "," : "]");
	}
}

static int guessrank(void) {
	static int init;
	static arb_t t;
	double ratio;

	if (!init) {
		arb_init(t);
		init = 1;
	}
	arb_div(t,acb_realref(bigdata[2]),acb_realref(bigdata[1]),prec);
	ratio = arf_get_d(arb_midref(t),ARF_RND_NEAR);
	if (epsilon < 0) {
		if (ratio <= 1) return 1;
		return 1+2*(int)floor(log(ratio/2)/log(4.0)+0.5);
	} else {
		if (ratio <= 1) return 0;
		return 2*(int)floor(log(ratio)/log(4.0)+0.5);
	}
}

// bisect a positive-to-negative sign change
static void bisect_zero(arb_t z,const arb_t t1in,arb_t t2in,int iterations) {
	static int init;
	static arb_t t,t1,t2,y;
	int i;

	if (!init) {
		arb_init(t);
		arb_init(t1);
		arb_init(t2);
		arb_init(y);
		init = 1;
	}
	arb_set(t1,t1in);
	arb_set(t2,t2in);
	for (i=0;i<iterations;i++) {
		arb_add(t,t1,t2,prec);
		arb_mul_2exp_si(t,t,-1);
		Lambda_adjusted(y,t);
		if (arb_contains_zero(y))
			break;
		if (arb_contains_positive(y))
			arb_swap(t,t1);
		else
			arb_swap(t,t2);
	}
	arb_add(t,t1,t2,prec);
	arb_mul_2exp_si(z,t,-1);
}

static int cmp_ui(const arb_t x,unsigned long int y) {
	static int init;
	static arb_t diff;

	if (!init) {
		arb_init(diff);
		init = 1;
	}
	arb_sub_ui(diff,x,y,prec);
	if (arb_is_negative(diff)) return -1;
	if (arb_is_positive(diff)) return 1;
	return 0;
}

static int cmp(const arb_t x,const arb_t y) {
	static int init;
	static arb_t diff;

	if (!init) {
		arb_init(diff);
		init = 1;
	}
	arb_sub(diff,x,y,prec);
	if (arb_is_negative(diff)) return -1;
	if (arb_is_positive(diff)) return 1;
	return 0;
}

static inline void avg(arb_t z,const arb_t x,const arb_t y) {
	arb_add(z,x,y,prec);
	arb_mul_2exp_si(z,z,-1);
}

static void bisect_min(arb_t z,arb_srcptr in[5],int iterations) {
	static int init;
	static arb_t t1,t2,t3,v1,v2,v3,x1,x2,y1,y2;
	int i;

	if (!init) {
		arb_init(t1); arb_init(t2); arb_init(t3);
		arb_init(v1); arb_init(v2); arb_init(v3);
		arb_init(x1); arb_init(x2);
		arb_init(y1); arb_init(y2);
		init = 1;
	}
	arb_set(t1,in[0]); arb_set(t3,in[1]); avg(t2,t1,t3);
	arb_set(v1,in[2]); arb_set(v2,in[3]); arb_set(v3,in[4]);
	for (i=0;i<iterations;i++) {
		avg(x1,t1,t2); avg(x2,t2,t3);
		Lambda_adjusted(y1,x1); Lambda_adjusted(y2,x2);
		if (cmp(v1,y1) > 0 && cmp(v2,y1) > 0) {
			arb_set(t3,t2); arb_set(t2,x1);
			arb_set(v3,v2); arb_set(v2,y1);
		} else if (cmp(v2,y2) > 0 && cmp(v3,y2) > 0) {
			arb_set(t1,t2); arb_set(t2,x2);
			arb_set(v1,v2); arb_set(v2,y2);
		} else {
			arb_set(t1,x1); arb_set(t3,x2);
			arb_set(v1,y1); arb_set(v3,y2);
		}
	}
	arb_set(z,t2);
}

static void upsample(void) {
	static int init;
	static arb_t t;
	int i,j;

	if (!init) {
		arb_init(t);
		init = 1;
	}

	for (i=0;i<sample_points;i++) {
		acb_set(bigdata[i*upsample_factor],data[i]);
		for (j=1;j<upsample_factor;j++)
			acb_zero(bigdata[i*upsample_factor+j]);
	}
	fft(bigdata);
	for (i=sample_points/2;i<big_sample_points-sample_points/2;i++)
		acb_zero(bigdata[i]);
	ifft(bigdata);
	for (i=0;i<big_sample_points;i++)
		acb_mul_ui(bigdata[i],bigdata[i],upsample_factor,prec);

	rank = guessrank();
	for (i=0;i<fpoints;i++) {
		arb_set_ui(t,i+1);
		arb_div_ui(t,t,bigA,prec);
		arb_pow_ui(t,t,rank,prec);
		// flip sign for ranks 2 and 3 so that leading order term is positive
		if (rank & 2) arb_neg(t,t);
		arb_div(f[i],acb_realref(bigdata[i+1]),t,prec);
	}
	nzeros = 0;
}

static void addzero(const arb_t z) {
	static int init;
	static arb_t t;
	int m;

	if (!init) {
		arb_init(t);
		init = 1;
	}

	if (nzeros == mzeros) {
		mzeros = mzeros ? mzeros*2 : 64;
		zero = (arb_t *)realloc(zero,mzeros*sizeof(arb_t));
		for (m=nzeros;m<mzeros;m++)
			arb_init(zero[m]);
	}
	arb_set(zero[nzeros++],z);

	for (m=0;m<fpoints;m++) {
		arb_set_ui(t,m+1);
		arb_div_ui(t,t,bigA,prec);
		arb_div(t,t,z,prec);
		arb_mul(t,t,t,prec);
		arb_neg(t,t); arb_add_ui(t,t,1,prec);
		arb_div(f[m],f[m],t,prec);
	}
}

static int findonezero(void) {
	static int init;
	static arb_t t1,t2,z;
	int m;
	arb_srcptr param[5];

	if (!init) {
		arb_init(t1);
		arb_init(t2);
		arb_init(z);
		init = 1;
	}

	// look for sign change
	for (m=0;m<fpoints;m++)
		if (arb_contains_zero(f[m])) {
			fprintf(stderr,"%ld: ambiguous sign\n",hash);
			exit(1);
		} else if (arb_is_negative(f[m])) {
			arb_set_ui(t1,m);
			arb_div_ui(t1,t1,bigA,prec);
			arb_set_ui(t2,m+1);
			arb_div_ui(t2,t2,bigA,prec);
			bisect_zero(z,t1,t2,prec);
			addzero(z);
			return 1;
		}

	// look for local minimum
	for (m=1;m<fpoints-1;m++)
		if (cmp(f[m-1],f[m]) > 0 && cmp(f[m+1],f[m]) > 0) {
			arb_set_ui(t1,m); arb_div_ui(t1,t1,bigA,prec);
			arb_set_ui(t2,m+2); arb_div_ui(t2,t2,bigA,prec);
			param[0] = t1; param[1] = t2;
			param[2] = f[m-1]; param[3] = f[m]; param[4] = f[m+1];
			bisect_min(t2,param,prec);
			Lambda_adjusted(z,t2);
			if (arb_is_negative(z)) {
				bisect_zero(z,t1,t2,prec);
				addzero(z);
				return 1;
			}
			arb_set_ui(t1,1);
			arb_mul_2exp_si(t1,t1,50-prec);
			if (cmp(z,t1) < 0) {
				addzero(t2); addzero(t2);
				return 2;
			} else {
				fprintf(stderr,"%ld: failed to isolate sign change\n",hash);
				exit(1);
			}
		}

	return 0;
}

int zero_cmp(const void *x,const void *y) {
	return cmp((arb_srcptr)x,(arb_srcptr)y);
}

void printzeros(void) {
	int m;

	while (findonezero());
	qsort(zero,nzeros,sizeof(zero[0]),zero_cmp);

	for (m=0;m<nzeros;m++)
		if (cmp_ui(zero[m],20) > 0)
			nzeros = m;

	printstr("[");
	for (m=nzeros-1;m>=0;m--) {
		printstr("-"); print(zero[m]); printstr(",");
	}
	for (m=rank;m>0;m--)
		printstr("0,");
	for (m=0;m<nzeros;m++) {
		print(zero[m]);
		printstr(m<nzeros-1 ? "," : "]");
	}
}

#define nprocs 64
int main(int argc,char *argv[]) {
	int i,k,d,line,kprocs=0;
	char buf[256],lpolys[256],*PQ,*s;
	long disc,cond,p;
	FILE *infile;
	arb_t c,t;

	if (argc != 2 || !(infile=fopen(argv[1],"r")) ) {
		printf("usage: %s infile > outfile\n",argv[0]);
		return 1;
	}

	data = fft_init(sample_points,prec);
	bigdata = fft_init(big_sample_points,prec);
	for (i=0;i<fpoints;i++) arb_init(f[i]);
	arb_init(c); arb_init(t); arb_init(scale);

	// do some initializations now so they
	// aren't repeated every time we fork
	nbuffer = 1<<20;
	buffer = (char *)calloc(nbuffer,1);
	arb_set_ui(t,1); K0(t,t);

	// scale = 2*Pi*A/sample_points
	arb_const_pi(scale,prec);
	arb_mul_ui(scale,scale,2*A,prec);
	arb_div_ui(scale,scale,sample_points,prec);

	line = 0;
	while (fgets(buf,sizeof(buf),infile)) {
		line++;
		if (!(s=strtok(buf,":")) || sscanf(s,"%ld",&disc) != 1) continue;
		if (!(s=strtok(NULL,":")) ) continue;
		if (strchr(s,'?') || sscanf(s,"%ld",&cond) != 1 || !cond)
			continue;
		if (!(s=strtok(NULL,":")) || sscanf(s,"%ld",&hash) != 1) continue;
		if (!(PQ=strtok(NULL,":"))) continue;
		if (!(s=strtok(NULL,":\n"))) continue;
		if (!strcmp(s,"1") || !strcmp(s,"-1")) {
			sscanf(s,"%d",&epsilon);
			if (!(s=strtok(NULL,":")) || !strcpy(lpolys,s)) continue;
			if (!(s=strtok(NULL,"\n"))) continue;
		} else
			continue;
		if ( !(curve=smalljac_curve_init(PQ,&k)) )
			continue;

		p = 2; k = 0;
		for (s=strtok(lpolys,",");s;s=strtok(NULL,","))
			for (;p<=disc;p++)
				if (disc % p == 0) {
					do disc /= p; while (disc % p == 0);
					bad_lfactors[k].p = p;
					d = i_poly_parse(bad_lfactors[k].f,4,s);
					for (i=d+1;i<=4;i++)
						bad_lfactors[k].f[i] = 0;
					k++;
					break;
				}
		bad_lfactors[k].p = 0;
		if (disc != 1) {
			fprintf(stderr,"something has gone wrong with bad L-factors\n");
			return 1;
		}

		if (!fork()) {
			compute_coeffs(sqrt((double)cond));
			arb_log_ui(c,cond,prec);
			arb_mul_2exp_si(c,c,-1);
			for (i=0;i<=sample_points/2;i++) {
				arb_mul_ui(t,scale,i,prec);
				arb_sub(t,c,t,prec);
				arb_exp(t,t,prec);
				S(t,t);

				// scale by 8 to get exactly complete L-function after fft
				arb_mul_2exp_si(t,t,3);
				acb_set_arb(data[i],t);
				if (i > 0 && i < sample_points/2)
					acb_mul_si(data[sample_points-i],data[i],epsilon,prec);
			}
			fft(data);
			if (epsilon < 0)
				for (i=0;i<sample_points;i++) {
					acb_mul_onei(data[i],data[i]);
					acb_neg(data[i],data[i]);
				}

			upsample();
			sprintf(buf,"%ld:%ld:%d:[[1,",hash,cond,epsilon);
			printstr(buf);
			if (rank > 0)
				printstr("0");
			else {
				arb_zero(t);
				gamma_factor(t,t);
				arb_div(t,scale,t,prec);
				arb_mul(t,t,acb_realref(data[0]),prec);
				print(t);
			}
			printstr("]]:");
			printzeros();
			printstr(":");
			printZ();
			printstr("\n");

			flock(1,LOCK_EX);
			write(1,buffer,strlen(buffer));
			flock(1,LOCK_UN);
			return 0;
		} else {
			if (kprocs == nprocs-1)
				wait(&i);
			else
				kprocs++;
		}
	}

	while (wait(&i) > 0);
	return 0;
}
