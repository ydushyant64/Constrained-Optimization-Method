#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int function_evaluations = 0; // Global variable to count function evaluations

double f(double x)
{
    function_evaluations++; // Increment the count each time f(x) is called
    return 20*sin(x) - 15*pow(x,2);
}

double df(double x) 
{
    double h = 1e-5;
    return (f(x + h) - f(x - h)) / (2 * h);
}

void bounding_phase(double *a, double *b, double x0, double delta) 
{
    int k = 0;
    int max_itr = 100;
    double x1;

    for (k = 1; k <= max_itr; k++) 
    {
        double f0 = f(x0);
        x1 = x0 - delta;  // Update x1
        double f1 = f(x1);

        if (f0 > f1) 
        {
            x0 = x1;  // Move x0 forward
            delta *= 2;  // Increase step size
        } 
        else 
        {
            *a = x0;
            *b = x1;
            printf("Number of iterations required in bounding phase: %d\n", k);
            printf("Maxima of the function lies between the region (%f, %f)\n", *a, *b);
            return;
        }
    }

    printf("Bounding phase did not converge within the maximum iterations.\n");
}

double secant_method(double a, double b)
 {
    double x1 = a;
    double x2 = b;
    double X1[100], X2[100];
    int i, max_itr = 100;
    double tol = 1e-5;
    double z;

    for (i = 0; i < max_itr; i++) 
    {
        z = x2 - df(x2) * (x2 - x1) / (df(x2) - df(x1));
        
        if (fabs(z - x2) < tol) 
        {
            printf("\nSecant method converged after %d iterations\n", i);
            break;
        }
        
        printf("\n%f\t%f",x1,x2);
        if (df(z) < 0) 
        {                         // Adjusting to ensure we approach maximum
            x2 = z;
        } 
        else 
        {
            x1 = z;
        }
        X1[i]=x1;
        X2[i]=x2;
    }

    printf("\nSecant method did not converge within the maximum iterations.\n");
    FILE *fp;
    fp = fopen("Result6.txt","w");
    for (int k = 0; k < i; k++) 
    {
    fprintf(fp,"%d\t%lf\t%lf\n",k,X1[k],X2[k]);
    }
    fclose(fp);
    return z;
}

int main() 
{
    double a=-4.0, b=4.0;
    int n =10;
    double delta = (b-a)/10.0;

int total_function_evaluations = 0;

    for (int i = 0; i < 10; i++) 
    {
        double x0;
        printf("Enter the initial guess between (-4,4) \n");
        scanf("%f",&x0);
        function_evaluations = 0; 
        bounding_phase(&a, &b, x0, delta);
        double min_x = secant_method(a, b);
        printf("Maxima found at x = %f with f(x) = %f\n", min_x, f(min_x));
        total_function_evaluations += function_evaluations;
    }

    printf("Total function evaluations over 10 runs: %d\n", total_function_evaluations);

    return 0;
}