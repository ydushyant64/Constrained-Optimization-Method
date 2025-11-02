#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include <time.h>


#define MAX_ITER 100
#define TOL 1e-6
#define epsilon 1e-3
#define epsilon2 1e-2

int counter =0;
int f_eval =0;

double objective_function(double x[]);
double pen_fun(double x[], int n, double R,double g[],int t,double P);
double single_var_pen_func(double alpha, double x[],double s[],int n, double R);
double gradient(double x[],double grad[], int n, double R, double g[],int t,double P); 
void Hassian(double x[], int n, double hassian[][n], double R, double g[],int t,double P);
double randInterval(double a1, double b1);
double bounding_phase(double x[], double s[], int n, double alpha,double R);
double secant_method(double x[], double lower_bound, double upper_bound, double s[], int n,double R);
double findBound(double x[], double S[], int n, double lower_alpha[], double upper_alpha[]);
double marquart(double x[], double lower_limit[], double upper_limit[], int n,int t, double R,double P,double g[]);
double mod(int n, double vector[]);




double pen_fun(double x[], int n, double R,double g[],int t,double P)
{
    counter++;
    f_eval++;

    //NOTE - Change input file name when changing objective function
    FILE *of = fopen("constraint violations.txt","a");

    // HB
    /*fprintf(of,"t\t\tR\t\t\t\tx1\t\t\t\tx2\t\t\t\tP\t\t\tg\n");
   
    g[0] = pow(x[0] - 5,2) + pow(x[1],2) - 26.0;
    
    if(g[0]>=0) // bracket operator returns only negative values
    {
        //printf("since g(0) is %lf\n",g[0]);
        g[0]=0;
        //printf("so g(0) is %lf\n",g[0]);
        
    }
    // if(fabs(g[0])<1e-2)
    // return 0;
    fprintf(of,"%d\t\t%0.6lf\t\t%0.6lf\t\t%0.6lf\t\t%0.6lf\t\t%0.6lf\n",t,R,x[0],x[1],P,g[0]);
    fclose(of);
    return pow(x[0]*x[0] + x[1] - 11,2) + pow(x[0] + x[1]*x[1] - 7,2) + R*(pow(g[0],2));
    */
     
    //1. Problem 1
    fprintf(of,"R \t\t x1 \t\t\t x2 \t\t\t g1 \t\t\t g2 \t\t\t\t P\n");
    double gnorm[2],g1max,g2max;
    g[0] = pow(x[0]-5,2) + pow(x[1]-5,2) - 100;
    g[1] = 82.81 - pow(x[0]-6,2) - pow(x[1]-5,2);
    if(g[0]>=0) // bracket operator returns only negative values
    {
        //printf("since g(0) is %lf\n",g[0]);
        g[0]=0;
        //printf("so g(0) is %lf\n",g[0]);
        
    }
    if(g[1]>=0)
    {
        //printf("since g(1) is %lf\n",g[1]);
        g[1]=0;
        //printf("so g(1) is %lf\n",g[1]);
    }
    //normalizing the constraints
    g1max = -35;
    g2max = -138.19;
    //printf("\nConstraint violation g1 = % lf \n",g1);
    //printf("Constraint violation g2 = % lf \n",g2);
    fprintf(of,"%0.1lf \t %lf \t\t %lf \t\t %lf \t\t %lf \t\t %lf\n",R,x[0],x[1],g[0],g[1],P);
    fclose(of);
    gnorm[0]=g[0]/g1max;
    gnorm[1]=g[1]/g2max;
    return pow(x[0]-10,3) + pow(x[1]-20,3) + R*(pow(gnorm[0],2) + pow(gnorm[1],2) );

    //2. Problem 2
    
    /*fprintf(of,"R \t\t x1 \t\t\t x2 \t\t\t g1 \t\t\t\t g2 \t\t\t P\n");
    double gnorm[2],g1max,g2max;
    g[0] = -1*pow(x[0],2) + x[1] - 1;
    g[1] = -1 + x[0] - pow(x[1]-4,2);
    if(g[0]>=0) // bracket operator returns only negative values
    {
        g[0]=0;
    }
    if(g[1]>=0)
    {
        g[1]=0;
    }
    //normalizing the constraints
    g1max = -101;
    g2max = -37;
    // printf("\nConstraint violation g1 = % lf \n",g1);
    // printf("Constraint violation g2 = % lf \n",g2);
    fprintf(of,"%0.1lf \t %lf \t\t %lf \t\t %lf \t\t %lf \t\t %lf\n",R,x[0],x[1],g[0],g[1],P);
    fclose(of);
    gnorm[0]=g[0]/g1max;
    gnorm[1]=g[1]/g2max;
    return -1*(pow( sin(2*3.142857*x[0]),3 )*sin(2*3.142857*x[1]))/( pow(x[0],3)*(x[0] + x[1]) ) + R*(pow(gnorm[0],2) + pow(gnorm[1],2) );
    */
    //3. Problem 3
    
    /*fprintf(of,"R \t\t g1 \t\t g2 \t\t g3 \t\t g4 \t\t g5 \t\t g6 \t\t\t\t x1 \t\t x2 \t\t\t\t x3 \t\t\t\t x4 \t\t\t\t x5 \t\t\t\t x6 \t\t\t\t x7 \t\t\t\t x8 \t\t\t\t P\n");
    double g1max,g2max,g3max,g4max,g5max,g6max,gnorm[6];
    g[0] = 1 - 0.0025*(x[3] + x[5]);
    g[1] = 1 - 0.0025*(-1*x[3] + x[4] + x[6]);
    g[2] = 1 - 0.01*(-1*x[5] + x[7]);
    g[3] = -100*x[0] + x[0]*x[5] - 833.33252*x[3] + 83333.333;
    g[4] = x[1]*(-1*x[3] + x[6]) + 1250*(x[3] - x[4]);
    g[5] = -1*x[2]*x[4] + x[2]*x[7] + 2500*x[4] - 1250000;
    for(int j=0;j<6;j++)
    {
        if(g[j]>=0) // bracket operator returns only negative values
        {
            g[j]=0;
        }
    }
    //normalizing the constraints
    g1max = -4;
    g2max = -3.975;
    g3max = -8.9;
    g4max = -1748999.187;
    g5max = -11137500;
    g6max = -11125000;
    fprintf(of,"%0.1lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",R,g[0],g[1],g[2],g[3],g[4],g[5],x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],P);
    fclose(of);
    gnorm[0]=g[0]/g1max;
    gnorm[1]=g[1]/g2max;
    gnorm[2]=g[2]/g3max;
    gnorm[3]=g[3]/g4max;
    gnorm[4]=g[4]/g5max;
    gnorm[5]=g[5]/g6max;
    return x[0] + x[1] + x[2] + R*(pow(gnorm[0],2) + pow(gnorm[1],2) + pow(gnorm[2],2) + pow(gnorm[3],2) + pow(gnorm[4],2) + pow(gnorm[5],2) );
    */
}

double single_var_pen_func(double alpha, double x[],double s[],int n, double R)
{
    //NOTE - Change input file name when changing objective function
    // HB
    /*double g1;
    g1 = pow(x[0]+alpha*s[0] - 5,2) + pow(x[1]+alpha*s[1],2) - 26.0;
    if(g1>=0)  // bracket operator returns only negative values
    {
        g1=0;
    }
    return (pow(x[0]+alpha*s[0] - 5,2) + pow(x[1]+alpha*s[1],2) - 26.0) +R*pow(g1,2);
*/

    //1. Problem 1
    double g1,g2,g1max,g2max;
    g1 = pow((x[0]+alpha*s[0])-5,2) + pow((x[1]+alpha*s[1])-5,2) - 100;
    g2 = 82.81 - pow((x[0]+alpha*s[0])-6,2) - pow((x[1]+alpha*s[1])-5,2);
    if(g1>=0)  // bracket operator returns only negative values
    {
        g1=0;
    }
    if(g2>=0)
    {
        g2=0;
    }
    g1max = -35;
    g2max = -138.19;
    g1=g1/g1max;
    g2=g2/g2max;

    return pow((x[0]+alpha*s[0])-10,3) + pow((x[1]+alpha*s[1])-20,3) + R*( pow(g1,2) + pow(g2,2) );

    //2. Problem 2
    
    /*double g1,g2,g1max,g2max;
    g1 = -1*pow((x[0]+alpha*s[0]),2) + (x[1]+alpha*s[1]) - 1;
    g2 = -1 + (x[0]+alpha*s[0]) - pow((x[1]+alpha*s[1])-4,2);
    if(g1>=0) // bracket operator returns only negative values
    {
        g1=0;
    }
    if(g2>=0)
    {
        g2=0;
    }
    //normalizing the constraints
    g1max = -101;
    g2max = -37;
    g1=g1/g1max;
    g2=g2/g2max;
    return -1*(pow( sin(2*3.142857*(x[0]+alpha*s[0])),3 )*sin(2*3.142857*(x[1]+alpha*s[1])))/( pow((x[0]+alpha*s[0]),3)*((x[0]+alpha*s[0]) + (x[1]+alpha*s[1])) ) + R*(pow(g1,2) + pow(g2,2) );
    */
    //3. Problem 3
    
    /*double g[6],g1max,g2max,g3max,g4max,g5max,g6max;
    g[0] = 1 - 0.0025*((x[3]+alpha*s[3]) + (x[5]+alpha*s[5]));
    g[1] = 1 - 0.0025*(-1*(x[3]+alpha*s[3]) + (x[4]+alpha*s[4]) + (x[6]+alpha*s[6]));
    g[2] = 1 - 0.01*(-1*(x[5]+alpha*s[5]) + (x[7]+alpha*s[7]));
    g[3] = -100*(x[0]+alpha*s[0]) + (x[0]+alpha*s[0])*(x[5]+alpha*s[5]) - 833.33252*(x[3]+alpha*s[3]) + 83333.333;
    g[4] = (x[1]+alpha*s[1])*(-1*(x[3]+alpha*s[3]) + (x[6]+alpha*s[6])) + 1250*((x[3]+alpha*s[3]) - (x[4]+alpha*s[4]));
    g[5] = -1*(x[2]+alpha*s[2])*(x[4]+alpha*s[4]) + (x[2]+alpha*s[2])*(x[7]+alpha*s[7]) + 2500*(x[4]+alpha*s[4]) - 1250000;
    for(int j=0;j<6;j++)
    {
        if(g[j]>=0) // bracket operator returns only negative values
        {
            g[j]=0;
        }
    }
    //normalizing the constraints
    g1max = -4;
    g2max = -3.975;
    g3max = -8.9;
    g4max = -1748999.187;
    g5max = -11137500;
    g6max = -11125000;
    g[0]=g[0]/g1max;
    g[1]=g[1]/g2max;
    g[2]=g[2]/g3max;
    g[3]=g[3]/g4max;
    g[4]=g[4]/g5max;
    g[5]=g[5]/g6max;
    return (x[0]+alpha*s[0]) + (x[1]+alpha*s[1]) + (x[2]+alpha*s[2]) + R*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2) + pow(g[3],2) + pow(g[4],2) + pow(g[5],2) );
    */
}

double objective_function(double x[])

{
    //HB
    //return pow(x[0]*x[0] + x[1] - 11,2) + pow(x[0] + x[1]*x[1] - 7,2);


    //1. Problem 1
    return pow(x[0]-10,3) + pow(x[1]-20,3) ;

    //2. Problem 2
    //return -1*(pow( sin(2*3.142857*x[0]),3 )*sin(2*3.142857*x[1]))/( pow(x[0],3)*(x[0] + x[1]) ) ;

    //3. Problem 3
    //return x[0] + x[1] + x[2] ;
}


    
double gradient(double x[],double grad[], int n, double R, double g[],int t,double P) 
{
    // Define a small value for numerical differentiation
    double h = 1e-3;
  
    // for(int i =0; i<n; i++)
    // {
    //     printf("Calculating gradient for x[%d] = %lf\n",i+1,x[i]);
    // }
    // Calculate partial derivatives using numerical differentiation
    for (int i = 0; i < n; i++) 
    {
        double x_temp = x[i];
        double f_temp = pen_fun(x,n,R,g,t,P);                        
        x[i] += h;
        grad[i] = (pen_fun(x,n,R,g,t,P) - f_temp) / h;
       // printf("Value of gradient at x[%d] is %lf\n",i+1,grad[i]);
        x[i] = x_temp;

    }
    return 0;
}

void Hassian(double x[], int n, double hassian[][n], double R, double g[],int t,double P)
{
    double h = 1e-3;
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++)     
        {
            double temp_i = x[i];
            double temp_j = x[j];
	    
            if (i==j)
            {                                              //will find diagonal elements
            	double f0 = pen_fun(x,n,R,g,t,P);
            	x[i] += h;
            	double f1 = pen_fun(x,n,R,g,t,P);

            	x[i] = temp_i;
            	x[i] -= h;
            	double f2 = pen_fun(x,n,R,g,t,P);
            	hassian[i][j] = ((f1 - 2*f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else
            {                                                   //will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = pen_fun(x,n,R,g,t,P);

            // Compute f(x - h, x_i)
            	x[i] = temp_i;
            	x[i] += h;
            	x[j] = temp_j;
            	x[j] -= h;
            	double f2 = pen_fun(x,n,R,g,t,P);

            // Compute f(x, x_j + h)
            	x[i] = temp_i;
            	x[i] -= h;
            	x[j] = temp_j;
            	x[j] += h;
            	double f3 = pen_fun(x,n,R,g,t,P);

            // Compute f(x, x_j - h)
            	x[i] = temp_i;
            	x[i] -=h;
            	x[j] = temp_j;
            	x[j] -=h;
            	double f4 = pen_fun(x,n,R,g,t,P);

            // Calculate the second partial derivative
            	hassian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	        }
            // Restore the original values
            x[i] = temp_i;
            x[j] = temp_j;
            //printf("\n ****************** Hassain Matrix *******************************************\n");
            printf("%10.3lf ", hassian[i][j]);
        }
        printf("\n");
    }
    //printf("\n");
}

double randNum()
{                                                             //will return random float between [0,1]
	return ( (double) rand() / (double) RAND_MAX);            //RAND_MAX is max random number generated by the machine
}

double randInterval(double a1, double b1)
{                 //will return random float between [a1,b1]
   double random;
	 random = randNum()*(b1-a1) + a1;
	return random;

}

double bounding_phase(double x[], double s[], int n, double alpha,double R)
{
    // a1 and a2 range of variables
    //x is a intitial point vector
    // x_new is a vector to store new point
    // s is vector to store the search direction

    double optimal;
      double lower_bound, upper_bound;
      double delta = 1e-2;
      double alpha_minus_delta = alpha - delta;
      double alpha_plus_delta = alpha + delta;

      
      int k =0;

      double f0, f0_minus_h, f0_plus_h;
      f0_minus_h =single_var_pen_func(alpha_minus_delta,x,s,n,R);
      f0 = single_var_pen_func(alpha,x,s,n,R);
      f0_plus_h = single_var_pen_func(alpha_plus_delta,x,s,n,R);

     if(f0_minus_h > f0 && f0 > f0_plus_h)
     {       // moving towards the right side (Minimization)
        while(f0 > f0_plus_h)
        {
            k++;
            alpha_minus_delta = alpha;
            alpha = alpha_plus_delta;
            alpha_plus_delta = alpha + pow(2,k)*delta;
            
            f0 = f0_minus_h;
            f0_minus_h = single_var_pen_func(alpha_minus_delta,x,s,n,R);
        }
        lower_bound = alpha_plus_delta;
        upper_bound = alpha_plus_delta;
     }
     else if(f0_minus_h < f0 && f0 < f0_plus_h)
     {    // moving towards the left side   (Maximization)
       
       while(f0 > f0_minus_h)
       {
        k++;
        alpha_plus_delta = alpha;
        alpha = alpha_minus_delta;
        alpha_minus_delta = alpha - pow(2,k)*delta;

        f0 = f0_minus_h;
        f0_minus_h = single_var_pen_func(alpha_minus_delta,x,s,n,R);
       }
       upper_bound = alpha_plus_delta;
       lower_bound = alpha_minus_delta;
     }
     else
     {
        for (int i =0; i<n; i++)
        {
            lower_bound = alpha_minus_delta;
            upper_bound = alpha_plus_delta;
        }
     }
     alpha = secant_method(x,lower_bound,upper_bound,s,n,R);
     //printf("Value of alpha is %lf\n",alpha);
     return alpha;
}
       


double secant_method(double x[], double lower_bound, double upper_bound, double s[], int n,double R)
{
 double h = 1e-3;
 double df2 = (single_var_pen_func(upper_bound + h, x, s, n,R) - single_var_pen_func(upper_bound - h, x, s, n,R)) / (2*h);
 double df1 = (single_var_pen_func(lower_bound + h, x, s, n,R) - single_var_pen_func(lower_bound - h, x, s, n,R)) / (2*h);  
 double new_alpha;
 for (int i = 0; i<MAX_ITER; i++)
 {
    new_alpha = upper_bound - (df2 *(upper_bound-lower_bound))/(df2 - df1);
    if(fabs(new_alpha - upper_bound) < TOL)
    {
       // printf("Secant Method converged after %d iterations\n",i);
        break;
    }
    double dz = (single_var_pen_func(new_alpha + h, x, s, n,R) - single_var_pen_func(new_alpha - h, x, s, n,R)) / (2*h);
    if(dz > 0)
    {
        upper_bound = new_alpha;
    }
    else
    {
        lower_bound = new_alpha;
    }
 }
 return upper_bound;
}

double findBound(double x[], double S[], int n, double lower_limit[], double upper_limit[])
{  
    double lower_alpha[n],upper_alpha[n];
    for(int i=0; i<n; i++)
    {
        lower_alpha[i] = lower_limit[i];
        upper_alpha[i] = upper_limit[i];
    } 
    double a[100],b[100];
    for(int i = 0; i < n; i++)
        {
            if((S[i]/mod(n,S)) > 0.0 )
            {
                a[i] = (lower_alpha[i] - x[i])*mod(n,S)/S[i];
                b[i] = (upper_alpha[i] - x[i])*mod(n,S)/S[i];
                
            }
            else
            {
                b[i] = (lower_alpha[i] - x[i])*mod(n,S)/S[i];
                a[i] = (upper_alpha[i] - x[i])*mod(n,S)/S[i];
            }
            
        }
    lower_alpha[0] = -100.0;
    upper_alpha[0] = 100.0;
    for(int i = 0; i < n; i++)
    {
        lower_alpha[0] = fmax(lower_alpha[i],a[i]);
        upper_alpha[0] = fmin(upper_alpha[i],b[i]);
    }
    return 0;
}

double mod(int n, double vector[])
{                         //vector of size n, will find its length/norm
	double sum =0;
	for(int i=0;i<n;i++)
    {
		sum = sum + pow(vector[i],2);
	}
	double modulus = sqrt(sum);
	return modulus;
}


void swapRows(int n, double A[][n], int row1, int row2) 
{
    for (int i = 0; i < n; i++) 
    {
        double temp = A[row1][i];     // Temporarily store the element at [row1][i]
        A[row1][i] = A[row2][i];      // Copy the element from [row2][i] to [row1][i]
        A[row2][i] = temp;            // Copy the original element from [row1][i] to [row2][i]
    }
}

void scaleRow(int n, double A[][n], int row, double factor) 
{
    for (int i = 0; i < n; i++) 
    {
        A[row][i] *= factor;          // Multiply each element in row 'row' by 'factor'
    }
}

void subtractRows(int n, double A[][n], int destRow, int srcRow, double factor) 
{
    for (int i = 0; i < n; i++) 
    {
        A[destRow][i] -= factor * A[srcRow][i];      // Subtract 'factor' * (srcRow element) from destRow element       
    }
}

int matrixInverse(int n, double A[][n], double B[n][n])     // A[][] is the hassian matrix + lamda*I
{
    double  C[n][n];                // B[][] will store the identity matrix    C[][] will store the copy of A[][] matrix

    // Copy A to C
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            C[i][j] = A[i][j];
        }
    }

    // creating an identity matrix B[][] 
   for(int i =0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (i==j)
            B[i][j] = 1;   // set the diagonal elements to 1
            else
            B[i][j] = 0;   // set the non diagonal elements to 0
        }
    }

    // Perform Gauss-Jordan elimination
    for (int i = 0; i < n; i++) 
    {
        // Find the pivot element and swap rows
        int pivotRow = i;
        while (C[pivotRow][i] == 0) 
        {
            pivotRow++;
            if (pivotRow == n) 
            {
                printf("Matrix is singular. Cannot find the inverse.\n");
                return 1; // Indicating error      
            }
        }
        swapRows(n, C, i, pivotRow);
        swapRows(n, B, i, pivotRow);

        // Scale the pivot row
        double pivot = C[i][i];
        scaleRow(n, C, i, 1.0 / pivot);
        scaleRow(n, B, i, 1.0 / pivot);

        // Eliminate other rows
        for (int j = 0; j < n; j++) 
        {
            if (j != i) 
            {
                double factor = C[j][i];
                subtractRows(n, C, j, i, factor);
                subtractRows(n, B, j, i, factor);
            }
        }
    }

    // Print the inverse matrix
    //printf("Inverse Matrix:\n");
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
           // printf("%lf\t", B[i][j]);
        }
        //printf("\n");
    }

   return 0; // Indicating success
}

double matrixMult(int n,double A[][n], double B[],double r[])
{          //will multiple a 2d array-A (size = n i.e. hessian inverse) with column vector - B (size =n, i.e. gradient) 
	for (int i=0;i<n;i++)
	{
		double sum=0;
		for (int j=0;j<n;j++)
		{
			sum = sum + A[i][j]*B[j];
		}
		r[i] = -sum;               // - negative sign added to satify the SEARCH direction requirement ( direction = -hessianInverse(A)*gradient(B))
	
	}
     //printf(" SEARCH DIRECTION \n");
    for(int i=0;i<n;i++)
    {
		//printf("Value at S[%d] is %.3lf\n",i,r[i]);
	}
}


double marquart(double x[], double lower_limit[], double upper_limit[], int n,int t, double R,double P,double g[])
{
    double I[n][n], lamda = 100.0, S[n], invH[n][n], Grad[n], alpha, x_1[n];
    double Norm_grad, x_diff[n];
    double INVERSE[n][n] , S_new[n];
    double A[n][n];   // A[n][n] stores the value of the hassian matrix H[n][n]
    double B[n][n];  // stores the inverse
    double C[n][n];  // copy of hassian
    // call the Hassian function
     
    double lower_alpha[n], upper_alpha[n];
    int iter;
    for(int i=0; i<n; i++)
    {
         lower_alpha[i] = lower_limit[i];
         upper_alpha[i] = upper_limit[i];
    }
    double func_val[100];
    printf("initial value of x are (%lf) , (%lf)\n",x[0],x[1]);
    // calculating S(0) to start the algorithm
    gradient(x,Grad,n,R,g,t,P);

    for(int i =0; i<n; i++)    // Identity Matrix 
    {
        for (int j=0; j<n; j++)
        {
            if (i==j)
            I[i][j] = 1;   // set the diagonal elements to 1
            else
            I[i][j] = 0;   // set the non diagonal elements to 0
        }
    }
    printf("\n --------------------------------- Hassain Matrix ------------------------------------------------\n");
    Hassian(x,n,A,R,g,t,P);
    // Add lamda to the hassan ( A + lamda*I)
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            {
                    invH[i][j] = A[i][j] + I[i][j]*lamda;
            }
    }
    matrixInverse(n,invH,INVERSE);
    matrixMult(n,INVERSE,Grad,S);

    double f_new, f_old =  pen_fun(x,n,R,g,t,P);
    printf("Old value of the function is %lf\n",f_old);
         
    // Unidirectional Search for alpha using the range of x
    findBound(x, S, n, lower_alpha, upper_alpha);
    alpha = lower_alpha[0] +  (upper_alpha[0] - lower_alpha[0])*(double) rand() / RAND_MAX ;
    double a_uni = lower_alpha[0], b_uni = upper_alpha[0];
    alpha = bounding_phase(x,S,n,alpha,R);
    alpha = 1.0;
            
    // calculation of x(1)
    for(int i = 0; i < n; i++)
    {
            x_1[i] = x[i] + alpha*S[i];
           // printf("x_new at x[%d] is %lf\n",i,x_1[i]);
    }

    for(int i = 0; i < n; i++)
    {
          x[i] = x_1[i];
    }
   
    for (iter =1; iter < MAX_ITER; iter++)
    {
        // Compute gradient of old and new x
        //gradient(x,Grad,n); 
        gradient(x,Grad,n,R,g,t,P);
        //printf("The Value of Gradient at iteraton %d is %0.3lf\n", iter,temp_grad);
        
        //Norm_grad = mod(n,Grad);
        //printf("Norm of the gradient is %lf\n",Norm_grad);
        printf("\n --------------------------------- Hassain Matrix ------------------------------------------------\n");
        Hassian(x,n,A,R,g,t,P);
        // Add lamda to the hassan ( A + lamda*I)
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                invH[i][j] = A[i][j] + I[i][j]*lamda;
            }
        }
        matrixInverse(n,invH,INVERSE);

        matrixMult(n,INVERSE,Grad,S);

        // //putting back values for next iteration
        //  for(int i = 0; i < n; i++)
        // {
        //   x[i] = x_1[i];
        //   S[i] = S_new[i];
        // }


        // f_old =  pen_fun(x,n,R,g,t,P);
        // printf("\nOld value of the function is %lf\n",f_old);
        // unidirectional search
        findBound(x, S, n, lower_alpha, upper_alpha);
        alpha = lower_alpha[0] +  (upper_alpha[0] - lower_alpha[0])*(double) rand() / RAND_MAX ;
        double a_uni = lower_alpha[0], b_uni = upper_alpha[0];
        alpha = bounding_phase(x,S,n,alpha,R);
        alpha = 1.0;
        // evaluating new x, for first iteration x(2)
        for(int i = 0; i < n; i++)
        {
            x_1[i] = x[i] + alpha*S[i];
            printf("x_new at x[%d] is %lf\n",i,x_1[i]);
        }

        f_new = pen_fun(x_1,n,R,g,t,P);
        printf("New value of the function is %lf\n",f_new);
        
        if(f_new <= f_old)
        {
            lamda = lamda/2;
            printf("Value of lambda is %lf\n",lamda);
        }
        else
        {
           //break;
            lamda = lamda*2;
            printf("Value of lamda is %lf\n",lamda);
            break;
            
        }
        
        
        //gradient(x_1,Grad,n);
        double G_norm = mod(n,Grad); //Evaluating Norm of new x, for first iteration

        // checking convergence
        if(G_norm < epsilon2)
        { 
            printf("\nFinal Value of x and x_new");
            
            for(int i = 0; i < n; i++)
            {   
                printf("\nx = %lf \t x_1 = %lf",x[i],x_1[i]);
                x[i] = x_1[i];
                
            }
            printf("\nNumber of itertions : %d\n",iter);
            FILE *fp1 = fopen("output_seq.txt","a");
            for(int i = 1 ; i < iter; i ++ )
                fprintf(fp1,"%d\t%lf\n",i,func_val[i]);  
            fclose(fp1);
            return f_new;
            break;
        }
        for(int i = 0; i < n; i++)
        {
          x[i] = x_1[i];
        }

        
        func_val[iter] = pen_fun(x,n,R,g,t,P);
        if(lamda==0.0)
        break;
    }
    FILE *fp1 = fopen("output_seq.txt","a");
    for(int i = 1 ; i < iter; i ++ )
    fprintf(fp1,"%d\t%lf\n",i,func_val[i]);  
    fclose(fp1);
            
    //fclose(fp1);
  printf("\nNumber of itertions : %d\n",iter);      
  return f_new;
}



int main()
{
int i,m,n,t,c;
FILE *table = fopen("Table.txt","a");
FILE *file1 = fopen("q1_input.txt","r");
    fscanf(file1,"no_of_variables %d\n",&n);
    fscanf(file1,"no_of_constraints %d\n",&m);
    double lower_limit[n],upper_limit[n],x[n],R,P_prev,P,g[m],x_initial[n];
 for(i=0;i<n;i++)
    {
        fscanf(file1,"lower_limit %*s %lf\n",&lower_limit[i]);
        fscanf(file1,"upper_limit %*s %lf\n",&upper_limit[i]);
        printf("lower limit for x(%d]) is %lf \tupper limit for x(%d) is %lf \n",i+1,lower_limit[i],i+1,upper_limit[i]);
    }
    fclose(file1);


//    n =2;
//    double x[2]={0.0,0.0},x_initial[2],R,P,P_prev,g[2];
//    double lower_limit[2]= {0.0, 0.0};
//    double upper_limit[2] = {5.0, 5.0};

srand(time(0));
//printf("Random variables for x1 and x2 are\n ");
printf("Rndom values of x1 x2 and x3 are\n");
for(i=0; i<n; i++)
{
x[i] = randInterval(lower_limit[i],upper_limit[i]);
printf("x[%d] = %lf\n",i,x[i]);
}
for(int i=0; i<n; i++)
{
    x_initial[i] = x[i];
}
 
// x[0] = 0.0;
// x[1] = 0.0;

printf("\n"); 

    t=0;
    c=10;
    R=0.1;
    // FILE *of = fopen("constraint violations.txt","a+");
    // fprintf(of,"R \t\t\t g1 \t\t\t g2\n");
    FILE *table2 = fopen("table2.txt","w");
     FILE *table3 = fopen("table3.txt","a");
    fprintf(table2,"sequence \t Function_evaluation\n");
    fprintf(table3,"Function Value \t Function Evaluation\n");
    do
    {
        if(t>0)
        {
            P_prev=P;
        }
        else
        {
            P_prev=0;
        }
        printf("\n\n****************************************************************************************************************************\n");
        printf("****************************************************************************************************************************\n");
        printf("\n Initializing Marquart's method taking R equal to %lf \n",R);
        printf("sequence no. %d\n\n",t);
        P=marquart(x,upper_limit,lower_limit,n,t,R,P,g);
        R=c*R;
        t++;
        // feval_total=feval_total+feval_seq;
        //fprintf(of,"%lf \t %lf \t %lf \n",R,g[0],g[1]);
        printf("\nfor sequence %d \t function evaluattion is %d\n",t-1,f_eval);
        fprintf(table2,"%d \t\t\t\t %d\n",t,f_eval);
        fprintf(table3,"%lf \t\t %d\n",P,f_eval);
        
        f_eval=0;
    }while(fabs(P-P_prev)>epsilon && R<pow(10,6));
    fclose(table2);
    fclose(table3);
   // fclose(of);
    // f=pen_fun(x,n,R,g);
    //fprintf(table,"Initial guess \t Optimal point \t sequnce \t f(x*) \t f_evaluation\n");
    
    for(int i=0; i<n; i++)
    {
    fprintf(table,"%lf, ",x_initial[i]);
    }
    for(int i=0; i<n; i++)
    {
    fprintf(table,"%lf, ",x[i]);
    }
    fprintf(table," \t %d \t %lf \t %d\n",t,P,counter);
    fclose(table);
    printf("\nAfter sequence number %d P(x) = %lf \n",t-1,P);\
    printf("\nTotal Function evaluations are %d\n",counter);
    return 0;
}


