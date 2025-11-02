#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include <time.h>


#define MAX_ITER 100
#define TOL 1e-6
#define epsilon 1e-6

int counter =0;

double objectiveFunction(double x[], double s[], double alpha, int n);
double single_variable(double x[], double alpha, double s[], int n);
double gradient(double x[],double grad[], int n);
void Hassian(double x[], int n, double hassian[][n]);
double randInterval(double a1, double b1);
double bounding_phase(double x[], double s[], int n, double alpha);
double secant_method(double x[],  double x1, double x2, double s[], int n);
double findBound(double x[], double S[], int n, double lower_alpha[], double upper_alpha[]);
double marquart(double x[], double lower, double upper, int n );
double mod(int n, double vector[]);


double objective_function(double x[], int n)
{
    counter ++;
    double sum = 0.0;
    //int i;

    // Himmelblau function
    
    //return pow((x[0]*x[0] + x[1] - 11),2) + pow((x[0] + x[1]*x[1] - 7),2);

    // // 1.Sum Squares Function          [-5.12 , 5.12]   , n = 5

    // for(int i=0; i<n; i++)
    // {
    //     sum = sum + (i+1)*pow(x[i],2.0);
    // }
    // return sum;

    // 2.Rosenbrock Function        [-2.048 , 2.048]  ,   n = 3

    for(int i=0; i<n-1; i++)
    {
        sum = sum + 100*pow((x[i+1]-pow(x[i],2.0)),2.0) + pow(x[i]-1,2.0);
    }
    return sum;

    // // 3.Dixon Price Function      [-10 , 10]     , n = 4    [0.707,1]

    // for(int i=1; i<n; i++)
    // {
    //     sum = sum + (i+1)*pow((2*pow(x[i],2.0)-x[i-1]),2.0);
    // }
    // return pow(x[0]-1,2.0) +sum;

    // // 4.Trid Function     [-n^2 , n^2]     , n = 6

    // for(int i=0; i<n; i++)
    // {
    //     sum = sum + pow(x[i]-1,2.0);
    // }
    // for(int i=1; i<n; i++)
    // {
    //     sum = sum - x[i]*x[i-1];
    // }
    // return sum;

    // // 5.Zakharov Function     [-5 , 10]   n =2

    // for(int i=0; i<n; i++)
    // {
    //     sum = sum + pow(x[i],2.0) + pow((0.5*(i+1)*x[i]),2) + pow((0.5*(i+1)*x[i]),4);
    // }  
    // return sum;
}


double single_variable(double x[], double alpha, double s[], int n)
{ //n = number of variables, alp =value of alpha ,x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION	
	double x_new[n];
	for (int i =0; i<n; i++)
	{            
		x_new[i] = x[i] + alpha*s[i];                
	}
	
	double f =  objective_function(x_new,n);
	return f;
}


    
double gradient(double x[],double grad[], int n) 
{
    // Define a small value for numerical differentiation
    double h = 1e-3;
  

    // Calculate partial derivatives using numerical differentiation
    for (int i = 0; i < n; i++) 
    {
        double x_temp = x[i];
        double f_temp = objective_function(x,n);
        x[i] += h;
        grad[i] = (objective_function(x,n) - f_temp) / h;
        printf("Value of gradient at x[%d] is %lf\n",i+1,grad[i]);
        x[i] = x_temp;

    }
    return 0;
}

void Hassian(double x[], int n, double hassian[][n])
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
            	double f0 = objective_function(x,n);
            	x[i] += h;
            	double f1 = objective_function(x,n);

            	x[i] = temp_i;
            	x[i] -= h;
            	double f2 = objective_function(x,n);
            	hassian[i][j] = ((f1 - 2*f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else
            {                                                   //will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = objective_function(x,n);

            // Compute f(x - h, x_i)
            	x[i] = temp_i;
            	x[i] += h;
            	x[j] = temp_j;
            	x[j] -= h;
            	double f2 = objective_function(x,n);

            // Compute f(x, x_j + h)
            	x[i] = temp_i;
            	x[i] -= h;
            	x[j] = temp_j;
            	x[j] += h;
            	double f3 = objective_function(x,n);

            // Compute f(x, x_j - h)
            	x[i] = temp_i;
            	x[i] -=h;
            	x[j] = temp_j;
            	x[j] -=h;
            	double f4 = objective_function(x,n);

            // Calculate the second partial derivative
            	hassian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	        }
            // Restore the original values
            x[i] = temp_i;
            x[j] = temp_j;
           // printf("\n ****************** Hassain Matrix *******************************************\n");
            printf("%10.3lf ", hassian[i][j]);
        }
        printf("\n");
    }
    printf("\n");
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

double bounding_phase(double x[], double s[], int n, double alpha)
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
      f0_minus_h = single_variable(x, alpha_minus_delta, s, n);
      f0 = single_variable(x, alpha, s, n);
      f0_plus_h = single_variable(x, alpha_plus_delta, s, n);

     if(f0_minus_h > f0 && f0 > f0_plus_h)
     {       // moving towards the right side (Minimization)
        while(f0 > f0_plus_h)
        {
            k++;
            alpha_minus_delta = alpha;
            alpha = alpha_plus_delta;
            alpha_plus_delta = alpha + pow(2,k)*delta;
            
            f0 = f0_minus_h;
            f0_minus_h = single_variable(x,alpha_minus_delta,s,n);
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
        f0_minus_h = single_variable(x,alpha_minus_delta,s,n);
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
     alpha = secant_method(x,lower_bound,upper_bound,s,n);
     printf("Value of alpha is %lf\n",alpha);
     return alpha;
}
       


double secant_method(double x[], double lower_bound, double upper_bound, double s[], int n)
{
 double h = 1e-3;
 double df2 = (single_variable(x, upper_bound + h, s, n) - single_variable(x,upper_bound - h, s, n)) / (2*h);
 double df1 = (single_variable(x, lower_bound + h, s, n) - single_variable(x, lower_bound - h, s, n)) / (2*h);  
 double new_alpha;
 for (int i = 0; i<MAX_ITER; i++)
 {
    new_alpha = upper_bound - (df2 *(upper_bound-lower_bound))/(df2 - df1);
    if(fabs(new_alpha - upper_bound) < TOL)
    {
        printf("Secant Method converged after %d iterations\n",i);
        break;
    }
    double dz = (single_variable(x, new_alpha + h, s, n) - single_variable(x, new_alpha -h, s, n)) / (2*h);
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

double findBound(double x[], double S[], int n, double lower_alpha[], double upper_alpha[])
{   
    double a[100],b[100];
    for(int i = 0; i < n; i++)
        {
            if((S[i]/mod(n,S)) > 0.0 )
            {
                a[i] = (lower_alpha[0] - x[i])*mod(n,S)/S[i];
                b[i] = (upper_alpha[0] - x[i])*mod(n,S)/S[i];
                
            }
            else
            {
                b[i] = (lower_alpha[0] - x[i])*mod(n,S)/S[i];
                a[i] = (upper_alpha[0] - x[i])*mod(n,S)/S[i];
            }
            
        }
    lower_alpha[0] = -100.0;
    upper_alpha[0] = 100.0;
    for(int i = 0; i < n; i++)
    {
        lower_alpha[0] = fmax(lower_alpha[0],a[i]);
        upper_alpha[0] = fmin(upper_alpha[0],b[i]);
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
    printf("Inverse Matrix:\n");
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            printf("%lf\t", B[i][j]);
        }
        printf("\n");
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
     printf(" SEARCH DIRECTION \n");
    for(int i=0;i<n;i++)
    {
		printf("Value at S[%d] is %.3lf\n",i,r[i]);
	}
}


double marquart(double x[], double lower, double upper, int n )
{
     double I[n][n], lamda = 100.0, S[n], invH[n][n], Grad[n], alpha, x_1[n];
     double Norm_grad, x_diff[n];
     double INVERSE[n][n] , S_new[n];
     double A[n][n];   // A[n][n] stores the value of the hassian matrix H[n][n]
     double B[n][n];  // stores the inverse
     double C[n][n];  // copy of hassian
     // call the Hassian function
     double lower_alpha[1], upper_alpha[1];
     int iter;
     lower_alpha[0] = lower;
     upper_alpha[0] = upper;
    double func_eval[100];
      // calculating S(0) to start the algorithm
      gradient(x,Grad,n);

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
       Hassian(x,n,A);
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

         double f_old =  objective_function(x,n);
         printf("Old value of the function is %lf\n",f_old);
         
         // Unidirectional Search for alpha using the range of x
        findBound(x, S, n, lower_alpha, upper_alpha);
        alpha = lower_alpha[0] +  (upper_alpha[0] - lower_alpha[0])*(double) rand() / RAND_MAX ;
        double a_uni = lower_alpha[0], b_uni = upper_alpha[0];
        alpha = bounding_phase(x,S,n,alpha);
        alpha = 1.0;
            
        // calculation of x(1)
        for(int i = 0; i < n; i++)
        {
            x_1[i] = x[i] + alpha*S[i];
            printf("x_new at x[%d] is %lf\n",i,x_1[i]);
        }

        for(int i = 0; i < n; i++)
        {
          x[i] = x_1[i];
        }
   
       for (iter =1; iter < MAX_ITER; iter++)
       {
        // Compute gradient of old and new x
        //gradient(x,Grad,n); 
        gradient(x,Grad,n);
        //printf("The Value of Gradient at iteraton %d is %0.3lf\n", iter,temp_grad);
        
        //Norm_grad = mod(n,Grad);
        //printf("Norm of the gradient is %lf\n",Norm_grad);
        Hassian(x,n,A);
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


        f_old =  objective_function(x,n);
        printf("Old value of the function is %lf\n",f_old);
        // unidirectional search
        findBound(x, S, n, lower_alpha, upper_alpha);
        alpha = lower_alpha[0] +  (upper_alpha[0] - lower_alpha[0])*(double) rand() / RAND_MAX ;
        double a_uni = lower_alpha[0], b_uni = upper_alpha[0];
        alpha = bounding_phase(x,S,n,alpha);
        alpha = 1.0;
        // evaluating new x, for first iteration x(2)
        for(int i = 0; i < n; i++)
        {
            x_1[i] = x[i] + alpha*S[i];
            printf("x_new at x[%d] is %lf\n",i,x_1[i]);
        }

        double f_new = objective_function(x_1,n);
        printf("New value of the function is %lf\n",f_new);
        
        if(f_new <= f_old)
        {
            lamda = lamda/2;
            printf("Value of lambda is %lf\n",lamda);
        }
        else
        {
            // //break;
            lamda = lamda*2;
            printf("Value of lamda is %lf\n",lamda);
            break;
            
        }
        
        //gradient(x_1,Grad,n);
        double G_norm = mod(n,Grad); //Evaluating Norm of new x, for first iteration

        // checking convergence
        if(G_norm < epsilon)
        { 
             printf("\nFinal Value of x and x_new");
            
            for(int i = 0; i < n; i++)
            {   
                printf("\nx = %lf \t x_1 = %lf",x[i],x_1[i]);
                x[i] = x_1[i];
                
            }
            break;
        }
         for(int i = 0; i < n; i++)
        {
          x[i] = x_1[i];
        }
        
        func_eval[iter] = objective_function(x,n);
  }
        FILE *fp1 = fopen("output_q1.txt","w");
        for(int i = 0 ; i < iter; i ++ )
            fprintf(fp1,"%d\t%lf\n",i+1,func_eval[i]);
            
    fclose(fp1);
        
  return iter;
}



int main()
{
double a, b;
int n;
//double x[2] = {0.0, 0.0};  
//Ask the user to input the number of random variables
printf("Enter the number of random variables: ");
scanf("%d", &n);
double x[n];
// Ask the user to input the lower and upper limits of the range
printf("Enter the lower limit of the range: ");
scanf("%lf", &a);
printf("Enter the upper limit of the range: ");
scanf("%lf", &b);


// // Seed the random number generator with the current time
srand(time(0));

// //Generate and print n random double values within the specified range
printf("Random variables: ");
for (int i = 0; i < n; i++) 
{
    x[i] = a +  (b - a)*(double) rand() / RAND_MAX ;
    printf("%lf\t",x[i]);  // Added space for readability
   
}
printf("\n");  // Added newline at the end
int iteration;
iteration  = marquart(x, a, b, n );

printf("Optimal Values are :\n ");
for(int i = 0; i < n; i++)
    printf("x[%d] : %lf\n",i+1,x[i]);

printf("\nOptimal function value : %lf\n",objective_function(x,n));

// FILE *fp1;
// fp1 = fopen("Table.txt","w");
// for (int i =0; i<10; i++)
// {
//     fprintf(fp1,"%d",i+1);
//     for(int j=0; j<n; j++)
//          fprintf(fp1,"\t%d",iter)
// }

printf("Total number of iterations %d\n",iteration);
    return 0;
}