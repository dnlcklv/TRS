// TRS2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <chrono>

using namespace std;

double F(double u)
{
    return u*u*u;
    //return 1;
}

double func(double x, double t)
{
    return ((1. - x) * (1. - x) - x * sin(t) - 2. * t);
}

double k1(double u)
{
    return sin(u);
    //return 1;
}

vector<double> TripleDiag(vector<double> A, vector<double> B, vector<double> C, vector<double> D) {
    vector<double> P(A.size());
    vector<double> Q(A.size());
    vector<double> X(A.size());
    int N = A.size();
    P[0] = C[0] / B[0];
    Q[0] = D[0] / B[0];
    for (int i = 1; i < N; ++i)
    {
        if (i < N - 1)
            P[i] = C[i] / (B[i] - A[i] * P[i - 1]);
        Q[i] = (D[i] - A[i] * Q[i - 1]) / (B[i] - A[i] * P[i - 1]);
    }

    // backward
    X[N - 1] = Q[N - 1];
    for (int i = N - 2; i >= 0; --i)
    {
        X[i] = Q[i] - P[i] * X[i + 1];
    }

    return X;
}

double u1_func(double x, double t)
{
    return x * cos(t) + t * (1 - x) * (1 - x);
}

void z1() {
    double ht, hx;

    vector<double> nX{ 10,100,130 };
    vector<double> nT(nX.size());

    for (int i = 0; i < nX.size(); i++) {
        clock_t start = clock();
        hx = 1. / (nX[i] - 1);
        nT[i] = 2 * (nX[i]) * (nX[i]);
        ht = 1. / (nT[i] - 1);
        vector<vector<double>> M(nT[i], vector<double>(nX[i])); // t , x ;
        ofstream file_to_cout("z1_" + to_string(i + 1) + ".txt");
        
        for (int k = 0; k < nX[i]; k++) {
            M[0][k] = hx * k;
        }
        double gamma = ht / hx / hx;
        for (int j = 1; j < nT[i]; j++) {

            for (int k = 0; k < nX[i]; k++) {
                if (k == 0) {
                    M[j][k] = j * ht;
                    continue;
                }
                if (k == nX[i] - 1) {
                    M[j][k] = (2. * cos(j * ht) * hx + M[j][k - 1]) / (1. + hx);
                    continue;
                }

                M[j][k] = gamma * (M[j - 1][k - 1] - 2. * M[j - 1][k] + M[j - 1][k + 1]) + ht * func(hx * k, ht * (j - 1)) + M[j - 1][k];
            }
        }

        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;

        for (int j = 0; j < nT[i]; j++) {
            for (int k = 0; k < nX[i]; k++) {
                file_to_cout << M[j][k] << " ";
            }
            file_to_cout << endl;
        }
        
        double dif;
        double max_dif;
        for (int j = 0; j < nT[i]; j++)
        {
            for (int k = 1; k < nX[i]; k++)
            {
                max_dif = 0;
                dif = abs(u1_func(k * hx, j * hx) - M[j][k]);
                if (dif > max_dif)
                {
                    max_dif = dif;
                }
            }
        }
        
        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
            << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size())  << "Time: " << time << endl;
    }
}

void z2_p1() {
    double ht, hx;
    vector<double> nX{ 10,100,130 };
    vector<double> nT = nX;
    for (int i = 0; i < nX.size(); i++) {

        clock_t start = clock();
        hx = 1. / (nX[i] - 1);
        //nT[i] = 2 * (nX[i]) * (nX[i]);
        ht = 1. / (nT[i] - 1);
        vector<vector<double>> M(nT[i], vector<double>(nX[i])); // t , x ;
        ofstream file_to_cout("z2_p1" + to_string(i + 1) + ".txt");
        vector<double> A(nX[i]), B(nX[i]), C(nX[i]), D(nX[i]);

        for (int k = 0; k < nX[i]; k++) {
            M[0][k] = k * hx;
        }
        A[0] = 0;
        B[0] = 1.;
        C[0] = 0;
        for (int k = 1; k < nX[i] - 1; k++) {
            A[k] = -1. / hx / hx;
            B[k] = 1. / ht + 2. / hx / hx;
            C[k] = -1. / hx / hx;
        }
        A[nX[i] - 1] = -1. / hx;
        B[nX[i] - 1] = 1. + 1. / hx;
        C[nX[i] - 1] = 0;
        for (int j = 1; j < nT[i]; j++) {
            D[0] = j * ht;
            for (int k = 1; k < nX[i] - 1; k++) {
                D[k] = func(hx * k, ht * j) + M[j - 1][k] / ht;
            }
            D[nX[i] - 1] = 2 * cos(ht * j);

            M[j] = TripleDiag(A, B, C, D);
        }

        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;

        for (int j = 0; j < nT[i]; j++) {//cout
            for (int k = 0; k < nX[i]; k++) {
                file_to_cout << M[j][k] << " ";
            }
            file_to_cout << endl;
        }

        double dif;
        double max_dif;
        for (int j = 0; j < nT[i]; j++)
        {
            for (int k = 1; k < nX[i]; k++)
            {
                max_dif = 0;
                dif = abs(u1_func(k * hx, j * hx) - M[j][k]);
                if (dif > max_dif)
                {
                    max_dif = dif;
                }
            }
        }

        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
            << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time: " << time << endl;
    }
}

void z2_p2() {
    double ht, hx;
    vector<double> nX{ 10,100,400 };
    vector<double> nT = nX;
    for (int i = 0; i < nX.size(); i++) {
        clock_t start = clock();
        hx = 1. / (nX[i] - 1);
        //nT[i] = 2 * (nX[i]) * (nX[i]);
        ht = 1. / (nT[i] - 1);
        vector<vector<double>> M(nT[i], vector<double>(nX[i])); // t , x ;
        ofstream file_to_cout("z2_p2" + to_string(i + 1) + ".txt");
        vector<double> A(nX[i]), B(nX[i]), C(nX[i]), D(nX[i]);

        for (int k = 0; k < nX[i]; k++) {
            M[0][k] = k * hx;
        }
        for (int j = 1; j < nT[i]; j++) {
            A[0] = 0;
            B[0] = 1;
            C[0] = 0;
            D[0] = j * ht;
            for (int k = 1; k < nX[i] - 1; k++) {
                A[k] = -1. / 2. / hx / hx;
                B[k] = 1. / ht + 1. / hx / hx;
                C[k] = -1. / 2. / hx / hx;
                D[k] = func(hx * k, ht * j) + 1. / 2. / hx / hx * (M[j - 1][k - 1] + M[j - 1][k + 1] - 2 * M[j - 1][k]) + 1. / ht * M[j - 1][k];
            }
            A[nX[i] - 1] = -1. / hx;
            B[nX[i] - 1] = 1 + 1. / hx;
            C[nX[i] - 1] = 0;
            D[nX[i] - 1] = 2 * cos(ht * j);
            M[j] = TripleDiag(A, B, C, D);
        }

        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;


        for (int j = 0; j < nT[i]; j++) {//cout
            for (int k = 0; k < nX[i]; k++) {
                file_to_cout << M[j][k] << " ";
            }
            file_to_cout << endl;
        }

        double dif;
        double max_dif;
        for (int j = 0; j < nT[i]; j++)
        {
            for (int k = 1; k < nX[i]; k++)
            {
                max_dif = 0;
                dif = abs(u1_func(k * hx, j * hx) - M[j][k]);
                if (dif > max_dif)
                {
                    max_dif = dif;
                }
            }
        }

        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
            << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time: " << time << endl;
    }
}

void z3()
{
    vector<double> nX {10,100,1000}; // по x
    vector<double> nT{ 10,100,1000};// по t
    for (int k = 0; k < nX.size(); k++)
    {
        ofstream result_file("z3_u" + to_string(k) + ".txt");
        auto begin = std::chrono::steady_clock::now();
        double ht = 1 / (nT[k] - 1);
        double hx = 1 / (nX[k] - 1);

        vector<vector<double>> u(nT[k], vector<double>(nX[k]));
        for (int i = 0; i < nX[k]; i++)
        {
            u[0][i] = i * hx; // началка
        }

        for (int j = 1; j < nT[k]; j++)
        {
            vector<double> A(nX[k]), B(nX[k]), C(nX[k]), D(nX[k]);
            //граничные
            A[0] = 0;
            B[0] = 1;
            C[0] = 0;
            D[0] = j * ht;
            A[nX[k] - 1] = -1. / hx;
            B[nX[k] - 1] = 1. + 1. / hx;
            C[nX[k] - 1] = 0;
            D[nX[k] - 1] = 2 * cos(j * ht);
            //
            double g = hx * hx / ht;
            for (int i = 1; i < nX[k] - 1; i++)
            {
                double ke = k1((u[j - 1][i + 1] + u[j - 1][i]) / 2);               
                double kw = k1((u[j - 1][i] + u[j - 1][i - 1]) / 2);
                A[i] = -kw / hx / hx;
                B[i] = 1. / ht + (kw + ke) / hx / hx;
                C[i] = -ke / hx / hx;
                D[i] = u[j - 1][i] / ht + func(i * hx,j * ht) * F(u[j - 1][i]);

            }
            u[j] = TripleDiag(A, B, C, D);
        }

        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

        for (int i = 0; i < nX[k]; i++)
        {
            for (int j = 0; j < nT[k]; j++)
            {
                result_file << u[j][i] << " ";
            }
            result_file << endl;
        }
        double dif;
        double max_dif;
        for (int i = 0; i < nT[k]; i++)
        {
            for (int j = 0; j < nX[k]; j++)
            {
                max_dif = 0;
                dif = abs(u1_func(j * hx, i * ht) - u[i][j]);
                if (dif > max_dif)
                {
                    max_dif = dif;
                }
            }
        }

        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
            << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time, ms: " << elapsed_ms.count() * 1e-3 << endl;
    }
} 

void z4(double eps)
{
    vector<double> nX{ 10,100 }; // по x
    vector<double> nT = nX;
    double delta = 0;
    for (int k = 0; k < nX.size(); k++)
    {
        ofstream result_file("z4_u" + to_string(k) + ".txt");

        auto begin = std::chrono::steady_clock::now();

        double ht = 1 / (nT[k] - 1);
        double hx = 1 / (nX[k] - 1);

        vector<vector<double>> u(nT[k], vector<double>(nX[k]));
        vector<double> previous;
        vector<double> A(nX[k]), B(nX[k]), C(nX[k]), D(nX[k]);

        for (int i = 0; i < nX[k]; i++)
        {
            u[0][i] = i * hx; // началка
        }

        for (int j = 1; j < nT[k]; j++)
        {
            //граничные
            A[0] = 0;
            B[0] = 1;
            C[0] = 0;
            D[0] = j * ht;
            A[nX[k] - 1] = -1. / hx;
            B[nX[k] - 1] = 1. + 1. / hx;
            C[nX[k] - 1] = 0;
            D[nX[k] - 1] = 2 * cos(j * ht);
            //
            double g = hx * hx / ht;
            for (int i = 1; i < nX[k] - 1; i++)
            {
                double ke = k1((u[j - 1][i + 1] + u[j - 1][i]) / 2);
                double kw = k1((u[j - 1][i] + u[j - 1][i - 1]) / 2);
                A[i] = -kw / hx / hx;
                B[i] = kw / hx / hx + ke / hx / hx + 1. / ht;
                C[i] = -ke / hx / hx;
                D[i] = func(i * hx, j * ht) * F(u[j - 1][i]) / hx + 1. / ht * u[j - 1][i];
            }
            u[j] = TripleDiag(A, B, C, D);



            int iter = 0;
            do
            {
                delta = 0;
                iter++;
                previous.assign(u[j].begin(), u[j].end());
                for (int i = 1; i < nX[k] - 1; i++)
                {
                    double ke = k1((u[j][i + 1] + u[j][i]) / 2);
                    double kw = k1((u[j][i] + u[j][i - 1]) / 2);
                    A[i] = -kw / hx / hx;
                    B[i] = kw / hx / hx + ke / hx / hx + 1. / ht;
                    C[i] = -ke / hx / hx;
                }
                u[j] = TripleDiag(A, B, C, D);

                for (int y = 1; y < nX[k]; y++)
                {
                    double diff = abs(k1(u[j][y] / 2 + u[j][y - 1] / 2) - k1(previous[y] / 2 + previous[y - 1] / 2));
                    if (diff > delta)
                    {
                        delta = diff;
                    }
                }
                previous.clear();
            } while (delta > eps);
            cout << "step x,t: " << hx << setw(13 + nX.size()) << "iter: " << iter << " error: " << delta << endl;
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

        for (int i = 0; i < nX[k]; i++)
        {
            for (int j = 0; j < nT[k]; j++)
            {
                result_file << u[j][i] << " ";
            }
            result_file << endl;
        }
       /* cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
            << setw(8 + nX.size()) << "Error: " << "..." << setw(8 + nX.size()) << "Time, ms: " << elapsed_ms.count() * 1e-3 << endl;*/

    }
}

int main()
{
    int ex;
    cout << "ex: " << endl;
    cin >> ex;
    switch (ex)
    {
    case 0:
    {
        z1();
        break;
    }
    case 1:
    {
        z2_p1();
        break;
    }
    case 2:
    {
        z2_p2();
        break;
    }
    case 3:
    {
        z3();
        break;
    }
    case 4: 
    {
        z4(1e-6);
        break;
    }
    }
}