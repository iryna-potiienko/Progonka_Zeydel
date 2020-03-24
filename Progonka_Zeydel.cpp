// Progonka_Zeydel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
using namespace std;
#include <Windows.h>
#include <cmath>

#define eps 0.001
#define n 3
double a[n][n] = { 3, 1, 0,
                  -1, 2, 0.5,
                   0, 0.5, 3 };

double b[n] = { 1, 1.75, 2.5 };

bool check() {
    double sum;
    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < n; j++) {
            if (i != j)sum += a[i][j];
        }
        if (a[i][i] < sum) return 1;
    }
    return 0;
}

double f(int i) {
    if (i < -1) return 0;
    if (i == -1) return 1;
    double det = a[i][i] * f( i - 1) - a[i + 1][i] * a[i][i + 1] * f( i - 2);
    return det;
}
double determinant() {
    int i;
    for (i = 0; i < n; i++) {
        f(i);
    }
    double det = f(i - 1);
    return det;
}
double matrix_norm(double b[][n]) {
    double sum, sumi[n] = {};
    for (int j = 0; j < n; j++) {
        sum = 0;
        for (int i = 0; i < n; i++) {
            sum += b[i][j];
        }
        sumi[j] = sum;
    }
    double max = 0;
    for (int j = 0; j < n; j++) {
        if (sumi[j] >= max) max = sumi[j];
    }
    return max;
}

double inversion()
{
    double A[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = a[i][j];
        }
    }
    double temp;

    double** E = new double* [n];

    for (int i = 0; i < n; i++)
        E[i] = new double[n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }
    }
    for (int k = 0; k < n; k++)
    {
        temp = A[k][k];

        for (int j = 0; j < n; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < n; i++)
        {
            temp = A[i][k];

            for (int j = 0; j < n; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = n - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];

            for (int j = 0; j < n; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = E[i][j];
        }
    }
    //printf("")
    for (int i = 0; i < n; i++) {
        printf("\n%f   %f   %f ", A[i][0], A[i][1], A[i][2]);
    }

    for (int i = 0; i < n; i++) {
        delete[] E[i];
    }

    delete[] E;
    double norm = matrix_norm(A);
    return norm;
}

double conol() {
    return inversion() * matrix_norm(a);
}

void progonka() {

    double alpha[n] = { }, beta[n] = {}, z[n];

    alpha[1] = -a[0][1]/a[0][0];
    beta[1] = b[0] / a[0][0];
    z[1] = a[1][1] + alpha[1] * a[1][0]; 

    alpha[2] = -a[1][2] / z[1];
    beta[2] = (b[1] - a[1][0] * beta[1]) / z[1];
    z[2] = a[2][2] + alpha[2] * a[2][1];

    for (int i = 1; i < n; i++) {
        cout << "\n alpha" << i << " = " << alpha[i] << " beta" << i << " = " << beta[i];
    }

    double y[n];
    y[n - 1] = (b[n - 1] - a[n - 1][n - 2] * beta[n - 1]) / (a[n - 1][n - 1] + alpha[n - 1] * a[n - 1][n - 2]);
    for (int i = n-1; i > 0; i--) {
        y[i - 1] = alpha[i] * y[i] + beta[i];
    } 
    
    cout << "\nВідповідь: ";
    for (int i = 0; i < n; i++) cout << "y" << i + 1 << " = " << y[i] << " ";
}

// Умова зупинки
bool converge(double xk[], double xkp[])
{
    double epsilon = 0.001;
    double norm = 0;
   // for (int i = 0; i < n; i++) cout << "\n " << xk[i]; //<<" "<< xkp[i];
    cout << "\n";
    for (int i = 0; i < n; i++)
        norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
    if (sqrt(norm) < eps) {
        return 0;
    }
    else return 1;
}

void zeidel() {
    double x[n], xp[n];
    if (check()!=0) {
        printf("\nНе виконується діагональна перевага!\n");
        return;
    }
    for (int j = 0; j < n; j++) x[j] = 0;
    do
    {
        for (int i = 0; i < n; i++)
            xp[i] = x[i];
        for (int i = 0; i < n; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (a[i][j] * x[j]);
            for (int j = i + 1; j < n; j++)
                var += (a[i][j] * xp[j]);
            x[i] = (b[i] - var) / a[i][i];
        }
        for (int i = 0; i < n; i++) cout << "\n x"<<i+1<<" = " << x[i];
    } while (converge(x, xp)!=0);
    
    cout << "Відповідь: ";
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << " ";
    }
}



int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
   
    printf("\nЗадана матриця:");
    for (int i = 0; i < n; i++) {
        printf("\n %f   %f   %f ", a[i][0], a[i][1], a[i][2]);
    }
    printf("\n\nМетод Зейделя: ");
    zeidel();
    printf("\n\nМетод прогонки: ");
    progonka();
    printf("\n\nВизначник: det(A) = %f \n", determinant());
    printf("\nОбернена матриця: ");
    //inversion(n);
    cout << "\n\nЧисло обумовленості: conol(A) = " << conol() << endl;
}


