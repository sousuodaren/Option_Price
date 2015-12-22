#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
int SIM_TIMES = 40000;
//#define SIM_TIMES		40000	
#define EXPIRATION_TIME         0.33
#define STRIKE_PRICE 	        100
#define INIT_PRICE		100

#define max(x,y) (x)>(y)?(x):(y)
#define CALL_OPTION 1 
#define A1 0.31938153
#define A2 -0.356563782
#define A3 1.781477937
#define A4 -1.821255978
#define A5 1.330274429
#define RSQRT2PI 0.39894228040143267793994605993438

pthread_mutex_t mutex; 
float S, K, T, sigma, r, R, mu, d;

static float * X;
float * Y1, * Y2, * Z1, * Z2;
int N, M;
int NE;
float dt, dh;
float th1=0.0, th2=0.0;
float c;
int Ps;

float *d_X;

float *d_Sum;

float d_Sum1,d_Sum2,d_Sum3,d_Sum4;
float *d_Y1,*d_Y2,*d_Z1,*d_Z2;

static float *Random_matrix;
static float *d_Random_matrix;
typedef struct argv_thread
{
	/* data */
	int i; //外层循环
	int a; //线程变化
};
int num_threads=10;

float MoroInvCND(float P){
    const float a1 = 2.50662823884f;
    const float a2 = -18.61500062529f;
    const float a3 = 41.39119773534f;
    const float a4 = -25.44106049637f;
    const float b1 = -8.4735109309f;
    const float b2 = 23.08336743743f;
    const float b3 = -21.06224101826f;
    const float b4 = 3.13082909833f;
    const float c1 = 0.337475482272615f;
    const float c2 = 0.976169019091719f;
    const float c3 = 0.160797971491821f;
    const float c4 = 2.76438810333863E-02f;
    const float c5 = 3.8405729373609E-03f;
    const float c6 = 3.951896511919E-04f;
    const float c7 = 3.21767881768E-05f;
    const float c8 = 2.888167364E-07f;
    const float c9 = 3.960315187E-07f;
    float y, z;

    if(P <= 0 || P >= 1.0f)
        return (float)(0x7FFFFFFF);

    y = P - 0.5f;
    if(fabsf(y) < 0.42f){
        z = y * y;
        z = y * (((a4 * z + a3) * z + a2) * z + a1) / ((((b4 * z + b3) * z + b2) * z + b1) * z + 1.0f);
    }else{
        if(y > 0)
            z = logf(-logf(1.0f - P));
        else
            z = logf(-logf(P));

        z = c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * (c6 + z * (c7 + z * (c8 + z * c9)))))));
        if(y < 0) z = -z;
    }

    return z;
}


float Ih(float y1, float y2, float x1, float x2, float x)
{
	float y;
	y = (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);
	return y;
}
float Ih1(float y1, float y2, float x1, float x2, float x)
{
	float y;
	y = (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);
	return y;
}
 float function_f(float y, float z, float mu, float sigma, float r, float d)
{
	float f;
	f = (-r) * y - 1 / sigma * (mu - r + d) * z;
	return f;
}
float function_f1(float y, float z, float mu, float sigma, float r, float d)
{
	float f;
	f = (-r) * y - 1 / sigma * (mu - r + d) * z;
	return f;
}
void Make_grid(float *X, int M, float dh)
{
	int i;
	for (i = 0; i <= M; i++)
	{
		X[i] = (i - M / 2) * dh;
	}
}

/*---Terminal_condition--*/
void Terminal_condition(int M, float * X, float * YT, float S0, float T, float K,
		float sigma, float mu, float r, float d)
{
	float St;
	int i;
        for (i = 0; i <= M; i++)
	{
		St = S0 * expf(sigma * X[i] + (mu - 0.5 * sigma * sigma) * T);

		if(CALL_OPTION == 1)
		{
			YT[i] = max(St-K,0.0);
		}
		else
		{	
			YT[i] = max(K-St,0.0);
		}
	}
}

void current_solution1(int j,  float *Y2,  float *Z2,float *Y1,  float *Z1,float *X,float th1, float th2, float dt, float dh, int NE, int N, float c, int M,float r, float sigma, float mu, float d,float *Random_matrix)
{
	int i, k, ii, a, Ps;
	
	float d_wt; 
	float Sy, Sz, Syw, Sf, Sfw; 
	float Xk; 
 	float Ey, Ez, Eyw, Ef, Efw;
	int size = (M + 1);
	unsigned int seed=1;
	float sq=sqrtf(dt);
	Ps = M / (2 * N);
	ii = Ps * (N - j); 

	for (i = ii; i <= M - ii; i++) 
	{
		
		Ey = Ez = Eyw = Ef = Efw = 0;
		for (k = 1; k <= NE; k++)
		{
			d_wt= Random_matrix[k];
			
			Xk = X[i] + d_wt;

			if (Xk < X[i - Ps]) 
				Xk = X[i - Ps]; 
			else if (Xk > X[i + Ps])
				Xk = X[i + Ps]; 

			a = (Xk - X[0]) / dh; 
			if (a == i + Ps) 
			{
				Sy = Y1[a];
				Sz = Z1[a];
			}
			else 
			{
				Sy = Ih1(Y1[a], Y1[a + 1], X[a], X[a + 1], Xk);
				Sz = Ih1(Z1[a], Z1[a + 1], X[a], X[a + 1], Xk);
			}

			Syw = Sy * d_wt;
			Sf = function_f1(Sy, Sz, mu, sigma, r, d);
			Sfw = Sf * d_wt;

			Ey = Sy + Ey;
			Ez = Sz + Ez;
			Eyw = Syw + Eyw;
			Ef = Sf + Ef;
			Efw = Sfw + Efw;
		}
		Z2[i] = (Eyw + dt * (1 - th2) * Efw - dt * (1 - th2) * Ez) / (NE * dt * th2);
		Y2[i] = ((Ey + dt * (1 - th1) * Ef) / NE- dt * th1 * (1 / sigma) * (mu - r + d) * Z2[i]) / (1 + dt * th1 * r);
	}
}
void sum(int i,int Ps,float *d_Y2,  float *d_Z2,float *d_Y1, float *d_Z1,float *d_X, float dh, int NE, float r, float sigma, float mu, float d,float *d_Random_matrix,float *d_Sum, int threadId)
{
		float Sy, Sz, Syw, Sf, Sfw; 
		float Ey, Ez, Eyw, Ef, Efw;
		float d_wt,Xk;
		int k,a;
		for (k = NE/num_threads*threadId; k < NE/num_threads*threadId+NE/num_threads; k++)
		{
		 	/* code */
		 	d_wt= d_Random_matrix[k];
			
			Xk = d_X[i] + d_wt;

			if (Xk < d_X[i - Ps]) 
				Xk = d_X[i - Ps]; 
			else if (Xk > d_X[i + Ps])
				Xk = d_X[i + Ps]; 

			a = (Xk - d_X[0]) / dh; 
			if (a == i + Ps) 
			{
				Sy = d_Y1[a];
				Sz = d_Z1[a];
			}
			else 
			{
				Sy = Ih(d_Y1[a], d_Y1[a + 1], d_X[a], d_X[a + 1], Xk);
				Sz = Ih(d_Z1[a], d_Z1[a + 1], d_X[a], d_X[a + 1], Xk);
			}

			Syw = Sy * d_wt;
			Sf = function_f(Sy, Sz, mu, sigma, r, d);
			Sfw = Sf * d_wt;

         //printf("%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", d_wt,Syw,Sfw,Sy,Ez);


			Ey = Sy + Ey;
			Ez = Sz + Ez;
			Eyw = Syw + Eyw;
			Ef = Sf + Ef;
			Efw = Sfw + Efw;
		 }
		 d_Sum[0]=Eyw+d_Sum[0];
		 d_Sum[1]=Efw+d_Sum[1];
		 d_Sum[2]=Ey+d_Sum[2];
		 d_Sum[3]=Ef+d_Sum[3];
}
void sumThread(void *argv)
{
	    struct argv_thread* i_a= (struct argv_thread*) argv; 
		int i=  i_a->i;
		int threadId = i_a->a;

		float Sy, Sz, Syw, Sf, Sfw; 
		Sy=Sz=Syw=Sf=Sfw=0;
		float Ey, Ez, Eyw, Ef, Efw;
		Ey=Ez=Eyw=Ef=Efw=0;
		float d_wt,Xk;
		int k,a;
		//printf("i is %d,threadId is %d\n",i,threadId );
		for (k = (NE/num_threads)*threadId; k <=(NE/num_threads)*threadId+(NE/num_threads); k++)
		{
		 	/* code */
		 	d_wt= d_Random_matrix[k];
			
			Xk = d_X[i] + d_wt;
//printf("k is %d,d_wt is %f,d_X[%d] is %f\n",k,d_wt,i,d_X[i] );
			if (Xk < d_X[i - Ps]) 
				Xk = d_X[i - Ps]; 
			else if (Xk > d_X[i + Ps])
				Xk = d_X[i + Ps]; 

//printf("%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", d_wt,Syw,Sfw,Sy,Ez);

			a = (Xk - d_X[0]) / dh; 
			if (a == i + Ps) 
			{
				Sy = d_Y1[a];
				Sz = d_Z1[a];
				//printf("******%10.4f,%10.4f\n",Sy,Ez);
			}else 
			{
				Sy = Ih(d_Y1[a], d_Y1[a + 1], d_X[a], d_X[a + 1], Xk);
				Sz = Ih(d_Z1[a], d_Z1[a + 1], d_X[a], d_X[a + 1], Xk);
			//	printf("++++++%10.4f,%10.4f\n",Sy,Ez);
			}

			Syw = Sy * d_wt;
			Sf = function_f(Sy, Sz, mu, sigma, r, d);
			Sfw = Sf * d_wt;

         //printf("%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", d_wt,Syw,Sfw,Sy,Ez);


			Ey = Sy + Ey;
			Ez = Sz + Ez;
			Eyw = Syw + Eyw;
			Ef = Sf + Ef;
			Efw = Sfw + Efw;
		 }
		 pthread_mutex_lock(&mutex);  
		 //printf("%d:%f,%f,%f,%f\n",threadId,Eyw,Efw,Ey,Ef);
		 d_Sum1=Eyw+d_Sum1;
		 d_Sum2=Efw+d_Sum2;
		 d_Sum3=Ey+d_Sum3;
		 d_Sum4=Ef+d_Sum4;
		 pthread_mutex_unlock(&mutex);  
}
//设备端调用，设备端执行
void current_solution(int j,  float *d_Y2,  float *d_Z2,float *d_Y1,  float *d_Z1,float *d_X,float th1, float th2, float dt, float dh, int NE, int N, float c, int M,float r, float sigma, float mu, float d,float *d_Random_matrix,float *d_Sum)
{
	int i, ii, Ps;
 
	int size = (M + 1);
	unsigned int seed=1;
	float sq=sqrtf(dt);
	Ps = M / (2 * N);
	ii = Ps * (N - j); 
	float Ey, Ez, Eyw, Ef, Efw;
	Ey = Ez = Eyw = Ef = Efw = 0;

	for (i = ii; i <= M - ii; i++) 
	{

	    d_Sum1=d_Sum2=d_Sum3=d_Sum4=0;
		//printf("%f,%f,%f,%f\n",d_Sum1,d_Sum2,d_Sum3,d_Sum4);
        
        int a;
		struct argv_thread argv[num_threads];
		pthread_mutex_init(&mutex, NULL);  
		pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
		for (a = 0; a < num_threads; a++) {
		    argv[a].i=i;
		    argv[a].a=a; 
			pthread_create(&pt[a], NULL, sumThread, (void *)&argv[a]);
		}
        for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);
        free(pt);
		//计算Eyw,Efw,Ey,Ef的累加和
		//for(threadId=0;threadId<200;threadId++)
		//	sum(i,Ps,d_Y2, d_Z2,d_Y1, d_Z1,d_X, dh,NE,r,sigma,mu,d,d_Random_matrix,d_Sum,threadId);
		//printf("%f,%f,%f,%f\n",d_Sum1,d_Sum2,d_Sum3,d_Sum4);
		Eyw=d_Sum1;
		Efw=d_Sum2;
		Ey=d_Sum3;
		Ef=d_Sum4;
//printf("%10.4f,%10.4f,%10.4f,%10.4f\n", Eyw,Efw,Ey,Ef);

		d_Z2[i] = (Eyw + dt * (1 - th2) * Efw - dt * (1 - th2) * Ez) / (NE * dt * th2);
		d_Y2[i] = ((Ey + dt * (1 - th1) * Ef) / NE- dt * th1 * (1 / sigma) * (mu - r + d) * d_Z2[i]) / (1 + dt * th1 * r);

	}
}
void submit(float *d_Y2,  float *d_Z2,float *d_Y1,  float *d_Z1,float *d_X,float th1, float th2, float dt, float dh, int NE, int N, float c, int M,float r, float sigma, float mu, float d,float *d_Random_matrix,float *d_Sum)
{
    int j;
	for (j = N - 1; j >= 0; j -= 2)
	{
		
		if (j == N - 1)
			th1 = th2 = 1;
		else
			th1 = th2 = 0.5;
			
		current_solution(j, d_Y2, d_Z2, d_Y1, d_Z1, d_X, th1, th2, dt, dh, NE, N, c, M, r,sigma, mu, d,d_Random_matrix,d_Sum);
		
		th1 = th2 = 0.5;
		
		if (j > 0)
		{
			current_solution(j - 1, d_Y1, d_Z1, d_Y2, d_Z2, d_X, th1, th2, dt, dh, NE, N,c, M, r, sigma, mu, d,d_Random_matrix,d_Sum);
		}
		else 
		{
			break;
		}
		printf("step.%d finish\n",j);
	}
}
void print_solution(float y, float z)
{
	printf("\n");
	printf("The value of Y0 is: %10.4f\n", y);
}
void print_array(float * A, int size)
{
	int i;
	for(i=0;i<size;i++){
		printf("A[%d] is %f ;", i,A[i]);	
	}
	printf("\n");
}


int main(int argc, char* argv[])
{
	

	int TIME_GRID = 64;

	if(argc>=2)
		TIME_GRID = atoi(argv[1]);
	if(argc>=3)
		SIM_TIMES = atoi(argv[2]);
	
	printf("TIME_GRID = %d\n",TIME_GRID);
	S = INIT_PRICE; 		
	K = STRIKE_PRICE; 	
	T = EXPIRATION_TIME; 

	sigma = 0.2;	
 	r = 0.03;
	R = 0.03;
	mu = 0.05;		
	d = 0.04;	
	
	N = TIME_GRID;
	dt = T / N;	

	dh = dt;	
	c = 5.0 * sqrtf(dt);	
	printf("c=%f\n",c);
	Ps = c / dh + 1;
	M = N * Ps * 2;	
	NE = SIM_TIMES;	


	int size = (M + 1) * sizeof(float);
	X = (float*) malloc(size);
//分配X的设备内存

	d_X=(float*) malloc(size);


	Y1 = (float*) malloc(size); 
	Y2 = (float*) malloc(size);
	Z1 = (float*) malloc(size); 
	Z2 = (float*) malloc(size);
//分配Y1,Y2,Z1,Z2的设备内存

	d_Y1=(float*) malloc(size);
	d_Y2=(float*) malloc(size);
	d_Z1=(float*) malloc(size);
	d_Z2=(float*) malloc(size);


	d_Sum=(float*) malloc(4*sizeof(float));


	memset(Y2,0,size);
	memset(Z2,0,size);
	memset(Z1,0,size);



	Make_grid(X, M, dh);
	
	Terminal_condition(M,  X,  Y1, S, T, K, sigma, mu, r, d);
	
	int j=0;
	
	struct timeval start;
	struct timeval finish;
	float tm;
	
	int k;
	unsigned int seed=(unsigned int )time(NULL);

	int num = NE+1;
	int rsize=sizeof(float)*(num);
	Random_matrix=(float *)malloc(rsize);

	d_Random_matrix=(float*) malloc(rsize);



	gettimeofday(&start,NULL);

    for(k=1;k<=num;k++)
    {
		Random_matrix[k]=MoroInvCND((float)k/(NE+1))*sqrt(dt);
	}
 //从主机端拷贝数据到设备端
	// memcpy(d_X,X,size);
	// memcpy(d_Y1,Y1,size);
	// memcpy(d_Y2,Y2,size);
	// memcpy(d_Z1,Z1,size);
	// memcpy(d_Z2,Z2,size);
	// memcpy(d_Random_matrix,Random_matrix,rsize);

//主机端调用，设备端执行
	submit(Y2,Z2,Y1,Z1,X, th1,  th2,  dt,  dh,  NE,  N,  c,  M, r,  sigma,  mu,  d, Random_matrix,d_Sum);
	
	gettimeofday(&finish,NULL); 
	tm=finish.tv_sec-start.tv_sec+(double)(finish.tv_usec-start.tv_usec)/1000000.0;
	printf("FINISHED!!\nALL Time is %.6f s\n",tm);
	if (j == -1){
		//memcpy(Y1,d_Y1,size);
		//memcpy(Z1,d_Z1,size);
		printf("%f,%f,%f,%f\n",Y1[M / 2], Z1[M / 2],d_Y1[M / 2], d_Z1[M / 2] );
		print_solution(d_Y1[M / 2], d_Z1[M / 2]);

	}
	else{
		//memcpy(Y2,d_Y2,size);
	   // memcpy(Z2,d_Z2,size);
	    printf("%f,%f,%f,%f\n",Y2[M / 2], Z2[M / 2],d_Y2[M / 2], d_Z2[M / 2] );
		print_solution(d_Y2[M / 2], d_Z2[M / 2]);

	}
		
	free(X);
	free(Random_matrix);
	free(Y1);
	free(Z1);
	free(Y2);
	free(Z2);

	free(d_X);
	free(d_Y1);
	free(d_Z1);
	free(d_Y2);
	free(d_Z2);
	free(d_Sum);
	free(d_Random_matrix);
}
