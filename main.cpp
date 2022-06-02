#include <iostream>
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

bool check(double **A)
{
    if (A[0][0] == 0)
    {
        cout << "Reference element is equal to 0 - Gauss method cannot be used!" << endl;
        return false;
    }
    return true;
}

void read(string filename, double **&A, double *&b, int &matrix_size)
{
    ifstream source_file(filename);
    if (!source_file.is_open())
    {
        cout << "The filename has not been open!" << endl;
    }
    source_file >> matrix_size;

    A = new double *[matrix_size];
    A[0] = new double[matrix_size * matrix_size];
    for (int i = 1; i < matrix_size; i++)
    {
        A[i] = A[i - 1] + matrix_size;
    }

    b = new double[matrix_size];
    char semicolumn;
    for (int i = 0; i < matrix_size + 1; i++)
    {
        source_file >> semicolumn;
    }

    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            source_file >> A[i][j];
            source_file >> semicolumn;
        }
        source_file >> semicolumn;
        source_file >> b[i];
    }
    source_file.close();
}

void showmatrix(double **A, double *b, int matrix_size)
{
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << b[i];
        cout << endl;
    }
}

void showresults(double *x, int matrix_size)
{
    for (int i = 0; i < matrix_size; i++)
    {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

void fillhistory(int *&history, int matrix_size)
{
    for (int i = 0; i < matrix_size; i++)
    {
        history[i] = i + 1;
    }
}

void fullselect(double **&A, double *&b, int matrix_size, int *&history, int startpoint)
{
    int mc = 0, mr = 0;
    double max = abs(A[startpoint][startpoint]);
    for (int i = startpoint; i < matrix_size; i++)
    {
        for (int j = startpoint; j < matrix_size; j++)
        {
            if (abs(A[i][j]) > max)
            {
                max = abs(A[i][j]);
                mr = i;
                mc = j;
            }
        }
    }
    if (mr != startpoint)
    {
        for (int i = startpoint; i < matrix_size; i++)
        {
            swap(A[startpoint][i], A[mr][i]);
        }
        swap(b[startpoint], b[mr]);
    }
    if (mc != startpoint)
    {
        for (int i = startpoint; i < matrix_size; i++)
        {
            swap(A[i][startpoint], A[i][mc]);
        }
        swap(history[startpoint], history[mc]);
    }
}

void rowselect(double **&A, double *&b, int *&history, int matrix_size, int startpoint)
{
    int mc = 0;
    double max = abs(A[startpoint][startpoint]);
    for (int i = startpoint; i < matrix_size; i++)
    {
        if (abs(A[startpoint][i]) > max)
        {
            max = abs(A[startpoint][i]);
            mc = i;
        }
    }
    if (mc != startpoint)
    {
        for (int i = startpoint; i < matrix_size; i++)
        {
            swap(A[i][startpoint], A[i][mc]);
        }
        swap(history[startpoint], history[mc]);
    }
}

void columnselect(double **&A, double *&b, int matrix_size, int startpoint)
{
    int mr = 0;
    double max = abs(A[startpoint][startpoint]);
    for (int i = 0; i < matrix_size; i++)
    {
        if (abs(A[i][startpoint]) > max)
        {
            max = abs(A[i][startpoint]);
            mr = i;
        }
    }
    if (mr != startpoint)
    {
        for (int i = startpoint; i < matrix_size; i++)
        {
            swap(A[startpoint][i], A[mr][i]);
        }
        swap(b[startpoint], b[mr]);
    }
}

void fixorder(double *&x, int *history, int matrix_size)
{
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size - 1; j++)
        {
            if (history[j] > history[j + 1])
            {
                swap(history[j], history[j + 1]);
                swap(x[j], x[j + 1]);
            }
        }
    }
}

void results(double **&A, double *&b, double *&x, int matrix_size)
{
    for (int i = matrix_size - 1; i >= 0; i--)
    {
        x[i] = b[i];
        for (int j = i + 1; j < matrix_size; j++)
        {
            x[i] -= A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
    }
}

void gauss(double **&A, double *&b, double *&x, int matrix_size)
{
    for (int i = 0; i < matrix_size - 1; i++)
    {
        for (int j = i + 1; j < matrix_size; j++)
        {
            double p = A[j][i] / A[i][i];
            for (int k = 0; k < matrix_size; k++)
            {
                A[j][k] = A[j][k] - p * A[i][k];
            }
            b[j] = b[j] - p * b[i];
        }
    }
    results(A, b, x, matrix_size);
}

void gaussfullselect(double **&A, double *&b, double *&x, int *&history, int matrix_size)
{
    for (int i = 0; i < matrix_size - 1; i++)
    {
        fullselect(A, b, matrix_size, history, i);
        for (int j = i + 1; j < matrix_size; j++)
        {
            double p = A[j][i] / A[i][i];
            for (int k = 0; k < matrix_size; k++)
            {
                A[j][k] = A[j][k] - p * A[i][k];
            }
            b[j] = b[j] - p * b[i];
        }
    }
    results(A, b, x, matrix_size);
}

void gaussrowselect(double **&A, double *&b, double *&x, int *&history, int matrix_size)
{
    for (int i = 0; i < matrix_size - 1; i++)
    {
        rowselect(A, b, history, matrix_size, i);
        for (int j = i + 1; j < matrix_size; j++)
        {
            double p = A[j][i] / A[i][i];
            for (int k = 0; k < matrix_size; k++)
            {
                A[j][k] = A[j][k] - p * A[i][k];
            }
            b[j] = b[j] - p * b[i];
        }
    }
    results(A, b, x, matrix_size);
}

void gausscolumnselect(double **&A, double *&b, double *&x, int matrix_size)
{
    for (int i = 0; i < matrix_size - 1; i++)
    {
        columnselect(A, b, matrix_size, i);
        for (int j = i + 1; j < matrix_size; j++)
        {
            double p = A[j][i] / A[i][i];
            for (int k = 0; k < matrix_size; k++)
            {
                A[j][k] = A[j][k] - p * A[i][k];
            }
            b[j] = b[j] - p * b[i];
        }
    }
    results(A, b, x, matrix_size);
}

void task1(string filename)
{
    cout << "==== 1 ====" << endl;
    double **A;
    double *b;
    double *x;
    int matrix_size;
    read(filename, A, b, matrix_size);
    x = new double[matrix_size];
    if (check(A))
    {
        gauss(A, b, x, matrix_size);
    }
    showresults(x, matrix_size);
    delete[] b;
    delete[] A[0];
    delete[] A;
    delete[] x;
}

void task2(string filename)
{
    cout << "==== 2 ====" << endl;
    double **A;
    double *b;
    double *x;
    int matrix_size;
    read(filename, A, b, matrix_size);
    x = new double[matrix_size];
    int *history;
    history = new int[matrix_size];
    fillhistory(history, matrix_size);
    gaussrowselect(A, b, x, history, matrix_size);
    fixorder(x, history, matrix_size);
    showresults(x, matrix_size);
    delete[] b;
    delete[] A[0];
    delete[] A;
    delete[] history;
    delete[] x;
}

void task3(string filename)
{
    cout << "==== 3 ====" << endl;
    double **A;
    double *b;
    double *x;
    int matrix_size;
    read(filename, A, b, matrix_size);
    x = new double[matrix_size];
    gausscolumnselect(A, b, x, matrix_size);
    showresults(x, matrix_size);
    delete[] b;
    delete[] A[0];
    delete[] A;
    delete[] x;
}

void task4(string filename)
{
    cout << "==== 4 ====" << endl;
    double **A;
    double *b;
    double *x;
    int matrix_size;
    read(filename, A, b, matrix_size);
    x = new double[matrix_size];
    int *history;
    history = new int[matrix_size];
    fillhistory(history, matrix_size);
    gaussfullselect(A, b, x, history, matrix_size);
    fixorder(x, history, matrix_size);
    showresults(x, matrix_size);
    delete[] b;
    delete[] A[0];
    delete[] A;
    delete[] history;
    delete[] x;
}

int main()
{
    task1("plik1.csv");
    task2("plik2.csv");
    task3("plik3.csv");
    task4("plik4.csv");
    return 0;
}
