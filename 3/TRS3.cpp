#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <Eigen/Sparse>,<Eigen/Dense>
#define pi 3.1415926535

using namespace std;
using namespace Eigen;

double u_analitic(double x, double y) 
{
	return sin(pi * x) * sin(pi * x) * sin(pi * y);
}

void createAnaliticFile(int N, int M)
{
	ofstream ftc("analitic.txt");
	double hx = 1. / (N - 1), hy = 1. / (M - 1);

	for (int i = 0; i <= N; i++) {
		ftc << u_analitic(hx * i, 0) << " ";//border y = 0;
	}
	ftc << endl;

	for (int j = 0; j < M - 1; j++) {
		for (int i = 0; i <= N; i++) {
			if (i == 0) {
				ftc << u_analitic(0, hy * j) << " ";//border x=0;
			}
			else if (i == N) {
				ftc << u_analitic(hx * i,hy * j) << " ";//border x = n;
			}
			else
				ftc << u_analitic(i * hx, j * hy) << " ";  // matrix_output
		}
		ftc << endl;
	}

	for (int i = 0; i <= N; i++) {
		ftc << u_analitic(hx * i, 1) << " "; //border y = n;
	}
}
double a(double x, double y)
{
	return exp(x);
}

double b(double x, double y)
{
	return exp(y);
}

double c(double x, double y)
{
	return 1;
}

double phi1(double y) { //x=0
	return 0;
}

double phi2(double y) { // x=lx
	return 0;
}

double phi3(double x) { // y=0
	return 0;
}

double phi4(double x) { // y=ly
	return 0;
}

double func(double x, double y) {
	return pi * pi * sin(pi * y) * (5 * cos(pi * x) * cos(pi * x) - 3);
}

double f_4(double x, double y) {
	return (-6.283185308) * exp(x) * sin(3.141592654 * x) * sin(3.141592654 * y) * cos(3.141592654 * x) - 19.73920881 * exp(x) * cos(3.141592654 * x) * cos(3.141592654 * x) * sin(3.141592654 * y) + 19.73920881 * exp(x) * sin(3.141592654 * x) * sin(3.141592654 * x) * sin(3.141592654 * y) - 3.141592654 * exp(y) * sin(3.141592654 * x) * sin(3.141592654 * x) * cos(3.141592654 * y) + 9.869604404 * exp(y) * sin(3.141592654 * x) * sin(3.141592654 * x) * sin(3.141592654 * y) - 1. * sin(3.141592654 * x) * sin(3.141592654 * x) * sin(3.141592654 * y);
}

double U_matrix_returning(int i, int j, int N, int M) {
	double hx = 1. / (N - 1);
	double hy = 1. / (M - 1);
	double alpha = -2. * (1. / hx / hx + 1. / hy / hy);
	double betta = 1. / hx / hx;
	double gamma = 1. / hy / hy;
	if (i == j) {
		return alpha;
	}
	else if (i - j == 1 && i % (N - 1) != 0) {
		return betta;
	}
	else if (i - j == -1 && j % (N - 1) != 0) {
		return betta;
	}
	else if (abs(i - j) == N - 1) {
		return gamma;
	}
	else
		return 0.;
}

void fill_nonlinear_matrix(int i, int j, int N, int M) {
	double hx = 1. / (N - 1);
	double hy = 1. / (M - 1);
	vector<double> left_gamma;
	vector<double> right_gamma;
	vector<double> left_betta;
	vector<double> right_betta;
	vector<double> alpha;

	for (int j = 1; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			left_gamma.push_back(-b(i * hx, j * hy) / 2. / hy + b(i * hx, j * hy) / hy / hy);
			right_gamma.push_back(b(i * hx, j * hy) / 2. / hy + b(i * hx, j * hy) / hy / hy);

			alpha.push_back(-2. * a(i * hx, j * hy) / hx / hx + b(i * hx, j * hy) / hy / hy + c(i * hx, j * hy));

			left_betta.push_back(-a(i * hx, j * hy) / 2 / hx + a(i * hx, j * hy) / hx / hx);
			right_betta.push_back(a(i * hx, j * hy) / 2 / hx + a(i * hx, j * hy) / hx / hx);
		}
	}
}

vector<double> Matrix_vector_Multi(int N, int M, vector<double>x) {
	vector<double> ax((N - 1) * (M - 1));
	double hx = 1. / (N - 1);
	double hy = 1. / (M - 1);
	double betta = 1. / hx / hx;
	double gamma = 1. / hy / hy;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			ax[j * (N - 1) + i] = 0;
			if (i > 0) {
				ax[j * (N - 1) + i] += betta * x[j * (N - 1) + i - 1]; // betta left
			}
			if (i < N - 2) {
				ax[j * (N - 1) + i] += betta * x[j * (N - 1) + i + 1]; //betta right
			}
			if (j < M - 2) {
				ax[j * (N - 1) + i] += gamma * x[(j + 1) * (N - 1) + i]; // gamma right
			}
			if (j > 0) {
				ax[j * (N - 1) + i] += gamma * x[(j - 1) * (N - 1) + i]; // gamma left
			}
			//*Idea: we are skippin all invalid memory operetions via including if-statements and summing all by parts:
		}
	}
	return ax;
}

void filling_F(vector<double>& f, int N, int M, int paramtoZadanie = 0) {
	double hx = 1. / (N - 1);
	double hy = 1. / (M - 1);
	if (paramtoZadanie != 0) {
		for (int j = 0; j < M - 1; j++) {
			for (int i = 0; i < N - 1; i++) {
				f[j * (N - 1) + i] = f_4((i + 1) * hx, (j + 1) * hy);
				if (j == 0) {
					f[j * (N - 1) + i] -= phi3((i + 1) * hx) / hy / hy;
				}
				if (j == (M - 1)) {
					f[j * (N - 1) + i] -= phi4((i + 1) * hx) / hy / hy;
				}
				if (i == (N - 1)) {
					f[j * (N - 1) + i] -= phi2((j + 1) * hy) / hx / hx;
				}
				if (i == 0) {
					f[j * (N - 1) + i] -= phi1((j + 1) * hy) / hx / hx;
				}
			}
		}
	}
	for (int j = 0; j < M - 1; j++) {
		for (int i = 0; i < N - 1; i++) {
			f[j * (N - 1) + i] = func((i + 1) * hx, (j + 1) * hy);
			if (j == 0) {
				f[j * (N - 1) + i] -= phi3((i + 1) * hx) / hy / hy;
			}
			if (j == (M - 1)) {
				f[j * (N - 1) + i] -= phi4((i + 1) * hx) / hy / hy;
			}
			if (i == (N - 1)) {
				f[j * (N - 1) + i] -= phi2((j + 1) * hy) / hx / hx;
			}
			if (i == 0) {
				f[j * (N - 1) + i] -= phi1((j + 1) * hy) / hx / hx;
			}
		}
	}

}
void filling_F(VectorXd& f, int N, int M, int paramtoZadanie = 0) {
	double hx = 1. / (N - 1);
	double hy = 1. / (M - 1);
	if (paramtoZadanie != 0) {
		for (int j = 0; j < M - 1; j++) {
			for (int i = 0; i < N - 1; i++) {
				f[j * (N - 1) + i] = -f_4((i + 1) * hx, (j + 1) * hy);
				if (j == 0) {
					f[j * (N - 1) + i] -= phi3((i + 1) * hx) / hy / hy;
				}
				if (j == (M - 1)) {
					f[j * (N - 1) + i] -= phi4((i + 1) * hx) / hy / hy;
				}
				if (i == (N - 1)) {
					f[j * (N - 1) + i] -= phi2((j + 1) * hy) / hx / hx;
				}
				if (i == 0) {
					f[j * (N - 1) + i] -= phi1((j + 1) * hy) / hx / hx;
				}
			}
		}
	}
	for (int j = 0; j < M - 1; j++) {
		for (int i = 0; i < N - 1; i++) {
			f[j * (N - 1) + i] = func((i + 1) * hx, (j + 1) * hy);
			if (j == 0) {
				f[j * (N - 1) + i] -= phi3((i + 1) * hx) / hy / hy;
			}
			if (j == (M - 1)) {
				f[j * (N - 1) + i] -= phi4((i + 1) * hx) / hy / hy;
			}
			if (i == (N - 1)) {
				f[j * (N - 1) + i] -= phi2((j + 1) * hy) / hx / hx;
			}
			if (i == 0) {
				f[j * (N - 1) + i] -= phi1((j + 1) * hy) / hx / hx;
			}
		}
	}

}

void LU_decompostion(int N, int M, vector<double> &L, vector<double> &U)
{
	int n = (N - 1) * (M - 1);

	for (int i = 0; i < n; i++)
	{
		L[i * n] = U_matrix_returning(i, 0, N, M);
		U[i] = U_matrix_returning(i, 0, N, M) / L[0];
	}

	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = U_matrix_returning(i, j,N,M);

			for (int k = 0; k < i; k++)
			{
				U[i * n + j] -= L[i * n + k] * U[k * n + j];
			}

			L[j * n + i] = U_matrix_returning(i, j, N, M);

			for (int k = 0; k < i; k++)
			{
				L[j * n + i] -= L[j * n + k] * U[k * n + i];
			}

			L[j * n + i] /= U[i * n + i];
		}
	}
}

void LU_decompostion1(int N, int M, vector<double>& L, vector<double>& U)
{
	int n = (N - 1) * (M - 1);

	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = U_matrix_returning(i, j, N, M);
		}
	}


	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			L[j * n + i] = U[j * n + i] / U[i * n + i];
		}
	}

	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++)
				L[j * n + i] = U[j * n + i] / U[i * n + i];

		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++)
				U[i * n + j] = U[i * n + j] - L[i * n + k - 1] * U[(k - 1) * n + j];
	}
}

vector<double> backward_up(int N, int M, vector<double> &U, vector<double> &b)
{
	int n = (N - 1) * (M - 1);

	vector<double> X(n * n);

	for (int i = n - 1; i >= 0; i--)
	{
		X[i] = b[i];
		for (int j = i + 1; j < n; j++)
		{
			X[i] -= U[i * n + j] * X[j];
		}
		X[i] /= U[i * n + i]; // /= должно быть
	}
	return X;
}

vector<double> backward_low(int N, int M, vector<double> &L, vector<double> &b)
{
	int n = (N - 1) * (M - 1);
	vector<double> Y(n * n);

	for (int i = 0; i < n; i++)
	{
		Y[i] = b[i];
		for (int j = 0; j < i; j++)
		{
			Y[i] -= L[i * n + j] * Y[j];
		}
		Y[i] /= L[i + i * n];
	}
	return Y;
}

vector<double> solve_LU(int N, int M, vector<double> &b)
{
	int n = (N - 1) * (M - 1);
	vector<double> L(n * n, 0), U(n * n, 0);
	LU_decompostion(N, M, L, U);
	auto Y = backward_low(N, M, L, b);
	auto X = backward_up(N, M, U, Y);

	return X;
}

double vec_vec_dot(vector<double>& rhs, vector<double>& lhs) {
	double result = 0;
	for (int i = 0; i < rhs.size(); i++) {
		result += rhs[i] * lhs[i];
	}
	return result;
}

vector<double> vectorSubtraction(vector<double> &A, vector<double> &B)
{
	vector<double> result(A.size(), 0);

	for (int i = 0; i < A.size(); i++)
	{
		result[i] = A[i] - B[i];
	}

	return result;
}

void z1_Poisson_problem_eq(int N, int M, int paramToCout = 0) {
	double hx = 1. / (N - 1), hy = 1. / (M - 1);
	vector<double> f((N - 1) * (M - 1));
	vector<double> x_prev((N - 1) * (M - 1), 0);
	double alpha = -2 * (1. / hx / hx + 1. / hy / hy);
	vector<double> x_next((N - 1) * (M - 1));
	for (int j = 0; j < M - 1; j++) {
		for (int i = 0; i < N - 1; i++) {
			f[j * (N - 1) + i] = func((i + 1) * hx, (j + 1) * hy);
			if (j == 0) {
				f[j * (N - 1) + i] -= phi3((i + 1) * hx) / hy / hy;
			}
			if (j == (M - 1)) {
				f[j * (N - 1) + i] -= phi4((i + 1) * hx) / hy / hy;
			}
			if (i == (N - 1)) {
				f[j * (N - 1) + i] -= phi2((j + 1) * hy) / hx / hx;
			}
			if (i == 0) {
				f[j * (N - 1) + i] -= phi1((j + 1) * hy) / hx / hx;
			}
		}
	}
	// x_k+1 = b - c*x_k; c = {aij,ii=0;}
	double diff = -1;
	x_next = f;

	vector<double> buf_Mult_MV;
	while (abs(diff) > 1e-15) {
		x_prev = x_next;
		diff = 0;
		for (int j = 0; j < (M - 1); j++) {
			buf_Mult_MV = Matrix_vector_Multi(N, M, x_prev);
			for (int i = 0; i < (N - 1); i++) {
				x_next[j * (N - 1) + i] = (f[j * (N - 1) + i] - buf_Mult_MV[j * (N - 1) + i]) / U_matrix_returning(i, i, N, M);
				diff = max(abs(x_prev[j * (N - 1) + i] - x_next[j * (N - 1) + i]), diff);
			}
		}
		//cout << diff << endl;
	}

	double diff_analitic = -1;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			diff_analitic = max(abs(u_analitic((i + 1) * hx, (j + 1) * hy) - x_next[j * (N - 1) + i]), diff_analitic); //err estimation
		}
	}
	cout << setprecision(3) << " hx: "  << hx << " hy: " << hy << "  Error: " << diff_analitic << endl;
	if (paramToCout == 1) {
		ofstream ftc("task1.txt");
		for (int i = 0; i <= N; i++) {
			ftc << phi3(hx * i) << " ";//border y = 0;
		}
		ftc << endl;

		for (int j = 0; j < M - 1; j++) {
			for (int i = 0; i <= N; i++) {
				if (i == 0) {
					ftc << phi1(hy * j) << " ";//border x=0;
				}
				else if (i == N) {
					ftc << phi2(hy * j) << " ";//border x = n;
				}
				else
					ftc << x_next[j * (N - 1) + i - 1] << " ";  // matrix_output
			}
			ftc << endl;
		}

		for (int i = 0; i <= N; i++) {
			ftc << phi4(hx * i) << " "; //border y = n;
		}
	}



}

void z2_Poisson_problem_eq(int N, int M, int paramToCout = 0) {
	double hx = 1. / (N - 1), hy = 1. / (M - 1);
	vector<double> f((N - 1) * (M - 1));
	vector<double> x_prev((N - 1) * (M - 1), 0);
	double alpha = -2 * (1. / hx / hx + 1. / hy / hy);
	double betta = 1. / hx / hx;
	double gamma = 1. / hy / hy;
	vector<double> x_next((N - 1) * (M - 1));

	filling_F(f, N, M);

	// x_k+1 = (1-omega)x_k  -  omega/aii (b - s1 -s2); c = {aij,ii=0;}
	double diff = -1;
	x_next = f;
	double omega = 1.5; //SOR parameter
	vector<double> buf_Mult_MV;
	while (abs(diff) > 1e-15) {
		x_prev = x_next;
		diff = 0;
		for (int j = 0; j < (M - 1); j++) {
			for (int i = 0; i < (N - 1); i++) {

				x_next[j * (N - 1) + i] = f[j * (N - 1) + i];
				if (i > 0) {
					x_next[j * (N - 1) + i] -= betta * x_next[j * (N - 1) + i - 1]; // betta left
				}
				if (i < N - 2) {
					x_next[j * (N - 1) + i] -= betta * x_prev[j * (N - 1) + i + 1]; //betta right
				}
				if (j < M - 2) {
					x_next[j * (N - 1) + i] -= gamma * x_prev[(j + 1) * (N - 1) + i]; // gamma right
				}
				if (j > 0) {
					x_next[j * (N - 1) + i] -= gamma * x_next[(j - 1) * (N - 1) + i]; // gamma left
				}
				x_next[j * (N - 1) + i] *= omega / U_matrix_returning(i, i, N, M);
				x_next[j * (N - 1) + i] += (1 - omega) * x_prev[j * (N - 1) + i];
				diff = max(abs(x_prev[j * (N - 1) + i] - x_next[j * (N - 1) + i]), diff);
			}
		}
		//cout << diff << endl;
	}

	double diff_analitic = -1;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			diff_analitic = max(abs(u_analitic((i + 1) * hx, (j + 1) * hy) - x_next[j * (N - 1) + i]), diff_analitic); //error estimation
		}
	}
	cout << setprecision(3) << " hx: " << hx << " hy: " << hy << "  Error: " << diff_analitic << endl;

	if (paramToCout == 1) {
		ofstream ftc("task2.txt");
		for (int i = 0; i <= N; i++) {
			ftc << phi3(hx * i) << " ";//border y = 0;
		}
		ftc << endl;

		for (int j = 0; j < M - 1; j++) {
			for (int i = 0; i <= N; i++) {
				if (i == 0) {
					ftc << phi1(hy * j) << " ";//border x=0;
				}
				else if (i == N) {
					ftc << phi2(hy * j) << " ";//border x = n;
				}
				else
					ftc << x_next[j * (N - 1) + i - 1] << " ";  // matrix_output
			}
			ftc << endl;
		}

		for (int i = 0; i <= N; i++) {
			ftc << phi4(hx * i) << " "; //border y = n;
		}
	}

}

void z3(int N, int M, int ptc)
{
	vector<double> f((N - 1) * (M - 1));
	double hx = 1. / (N - 1), hy = 1. / (M - 1);
	filling_F(f, N, M);
	auto solve = solve_LU(N, M, f);

	double diff_analitic = -1;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			diff_analitic = max(abs(u_analitic((i + 1) * hx, (j + 1) * hy) - solve[j * (N - 1) + i]), diff_analitic);
		}
	}
	cout << setprecision(3) << " hx: " << hx << " hy: " << hy << "  Error: " << diff_analitic << endl;

	if (ptc == 1) 
	{
		ofstream ftc("task3.txt");
		for (int i = 0; i <= N; i++) {
			ftc << phi3(hx * i) << " ";//border y = 0;
		}
		ftc << endl;

		for (int j = 0; j < M - 1; j++) {
			for (int i = 0; i <= N; i++) {
				if (i == 0) {
					ftc << phi1(hy * j) << " ";//border x=0;
				}
				else if (i == N) {
					ftc << phi2(hy * j) << " ";//border x = n;
				}
				else
					ftc << solve[j * (N - 1) + i - 1] << " ";  // matrix_output
			}
			ftc << endl;
		}

		for (int i = 0; i <= N; i++) {
			ftc << phi4(hx * i) << " "; //border y = n;
		}
	}
}

void z4_Poisson_problem_eq_notwork(int N, int M, int paramToCout = 0) {
	//cgm
	double hx = 1. / (N - 1), hy = 1. / (M - 1);
	int iter = 0;
	double alpha = 0;
	double betta = 0;
	double err = 0; double x_buf = 0;


	vector<double> r_prev((N - 1) * (M - 1));
	vector<double> f((N - 1) * (M - 1));
	filling_F(f, N, M, 1);
	filling_F(r_prev, N, M, 1);
	vector<double> x((N - 1) * (M - 1), 0.);
	vector<double> r_next((N - 1) * (M - 1));
	vector<double>Az((N - 1) * (M - 1));





	vector<double>z_prev = r_prev;



	double x_coord = 0;
	double y_coord = 0;
	do {
		err = 0;

		for (int j = 0; j < M - 1; j++) {
			for (int i = 0; i < N - 1; i++) {
				x_coord = hx * (i + 1);
				y_coord = hy * (j + 1);


				Az[j * (N - 1) + i] = (1 - 2 * (exp(x_coord) / hx / hx + exp(y_coord) / hy / hy)) * z_prev[j * (N - 1) + i];
				if (i > 0) {
					Az[j * (N - 1) + i] += exp(x_coord) * (1. / hx / hx - 1. / 2. / hx) * z_prev[j * (N - 1) + i - 1]; // betta left
				}
				if (i < N - 2) {
					Az[j * (N - 1) + i] += exp(x_coord) * (1. / hx / hx + 1. / 2. / hx) * z_prev[j * (N - 1) + i + 1]; //betta right
				}
				if (j < M - 2) {
					Az[j * (N - 1) + i] += exp(y_coord) * (1. / hy / hy + 1. / 2. / hy) * z_prev[(j + 1) * (N - 1) + i]; // gamma right
				}
				if (j > 0) {
					Az[j * (N - 1) + i] += exp(y_coord) * (1. / hy / hy - 1. / 2. / hy) * z_prev[(j - 1) * (N - 1) + i]; // gamma left
				}
			}
		}
		alpha = vec_vec_dot(r_prev, r_prev) / vec_vec_dot(Az, z_prev);
		for (int i = 0; i < x.size(); i++) {
			x_buf = x[i];
			x[i] += alpha * z_prev[i];
			r_next[i] = r_prev[i] - alpha * Az[i];
			err = max(err, abs(x[i] - x_buf));
		}
		betta = vec_vec_dot(r_next, r_next) / vec_vec_dot(r_prev, r_prev);
		for (int i = 0; i < x.size(); i++) {
			z_prev[i] = r_next[i] + betta * z_prev[i];
		}
		iter++;
		//cout << sqrt(vec_vec_dot(r_next, r_next)) / sqrt(vec_vec_dot(f, f)) << "\tIter:" << iter <<endl;
		r_prev = r_next;
	} while (err > 1e-15);

	double diff_analitic = -1;
	double max_num = 0;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			diff_analitic = max(abs(u_analitic((i + 1) * hx, (j + 1) * hy) - x[j * (N - 1) + i]), diff_analitic); //error estimation
			max_num = max(abs(x[j * (N - 1) + i]), max_num);
		}
	}
	//cout << "MCGM" << endl;
	cout << setprecision(3) << " hx: " << hx << " hy: " << hy << "  Error: " << diff_analitic << endl;
}


void z4_Eigen(int N, int M, int paramToCout = 0) {
	double hx = 1. / (N - 1), hy = 1. / (M - 1);
	VectorXd f((N - 1) * (M - 1));
	filling_F(f, N, M, 1);
	double value = 0;

	double x_coord, y_coord;
	vector<Triplet<double>> triplets;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			x_coord = hx * (i + 1);
			y_coord = hy * (j + 1);
			triplets.push_back(Triplet<double>(j * (N - 1) + i, j * (N - 1) + i, (1 - 2 * (exp(x_coord) / hx / hx + exp(y_coord) / hy / hy)))); // alpha
		}
	}
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 1; i < (N - 1); i++) {
			x_coord = hx * (i + 1);
			y_coord = hy * (j + 1);
			triplets.push_back(Triplet<double>(j * (N - 1) + i, j * (N - 1) + i - 1, exp(x_coord) * (1. / hx / hx - 1. / 2. / hx))); // betta left
		}
	}
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 2); i++) {
			x_coord = hx * (i + 1);
			y_coord = hy * (j + 1);
			triplets.push_back(Triplet<double>(j * (N - 1) + i, j * (N - 1) + i + 1, exp(x_coord) * (1. / hx / hx + 1. / 2. / hx))); // betta right
		}
	}
	for (int j = 0; j < (M - 2); j++) {
		for (int i = 0; i < (N - 1); i++) {
			x_coord = hx * (i + 1);
			y_coord = hy * (j + 1);
			triplets.push_back(Triplet<double>(j * (N - 1) + i, (j + 1) * (N - 1) + i, exp(y_coord) * (1. / hy / hy + 1. / 2. / hy))); // gamma right
		}
	}
	for (int j = 1; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			x_coord = hx * (i + 1);
			y_coord = hy * (j + 1);
			triplets.push_back(Triplet<double>(j * (N - 1) + i, (j - 1) * (N - 1) + i, exp(y_coord) * (1. / hy / hy - 1. / 2. / hy))); // gamma left
		}
	}

	SparseMatrix<double> A((N - 1) * (M - 1), (N - 1) * (M - 1));
	A.setFromTriplets(triplets.begin(), triplets.end());
	ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
	cg.compute(A);
	auto x = cg.solve(f);


	double diff_analitic = -1;
	double max_num = 0;
	for (int j = 0; j < (M - 1); j++) {
		for (int i = 0; i < (N - 1); i++) {
			diff_analitic = max(abs(u_analitic((i + 1) * hx, (j + 1) * hy) - x[j * (N - 1) + i]), diff_analitic); //error estimation
			max_num = max(abs(x[j * (N - 1) + i]), max_num);
		}
	}
	//cout << "CG" << endl;
	cout << setprecision(3) << " hx: " << hx << " hy: " << hy << "  Error: " << diff_analitic << endl;



}

int main() {
	int Number_of_ex;

	cin >> Number_of_ex;
	vector <int> N = { 4,10,20,30,40,50,60 };
	switch (Number_of_ex)
	{
	case 1:
		for (auto& el : N)
		{
			z1_Poisson_problem_eq(el, el);
		}
		//z1_Poisson_problem_eq(30, 30,1);
		break;
	case 2:
		for (auto& el : N)
		{
			z2_Poisson_problem_eq(el, el);
		}
		//z2_Poisson_problem_eq(30, 30, 1);
		break;
	case 3: 
		for (auto& el : N)
		{
			z3(el, el, 1);
		}
	case 4:
		for (auto& el : N)
		{
			//z4_Eigen(el, el);
			z4_Poisson_problem_eq_notwork(el, el);

			cout << "=====";
			//Jacobi_solver(el, el, 1);
		}
		//Jacobi_solver(30, 30, 1);
	case 5:
		createAnaliticFile(N[N.size() - 1], N[N.size() - 1]);
		break;
	}
}