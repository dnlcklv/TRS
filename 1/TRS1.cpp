#include <iostream>
#include <vector>
#include<string>
#include <fstream>

using namespace std;

double U(double x) {
	return -(4 / 9.1) * (x * x * sinh(x) + 2 * x * cosh(x));
}

double func(double x) {
	return 6 * x * x + 1;
};

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

double f1(double x, double v, double u)
{
	return 4 * v + 6 * x * x + 1;
}

double f2(double x, double v, double u)
{
	return v;
}

double rk4(double p)
{
	int n = 1000;
	double k1, k2, k3, k4;
	double q1, q2, q3, q4;
	double h = 1. / n;
	double v_0, u_n, v_n, u_0;
	u_0 = v_0 = 0;
	u_n = p;
	v_n = u_n - 3.5625;

	double x = 0;
	for (int i = 0; i < n; i++)
	{
		u_0 = u_n;
		v_0 = v_n;
		x += i * h;
		k1 = f1(x, v_0, u_0);
		q1 = f2(x, v_0, u_0);
		k2 = f1(x + h / 3, v_0 + h * k1 / 3, u_0 + h * q1 / 3);
		q2 = f2(x + h / 3, v_0 + h * k1 / 3, u_0 + h * q1 / 3);
		k3 = f1(x + 2 * h / 3, v_0 - h * k1 / 3, u_0 - h * q1 / 3 + h * q2);
		q3 = f2(x + 2 * h / 3, v_0 - h * k1 / 3, u_0 - h * q1 / 3 + h * q2);
		k4 = f1(x + h, v_0 + h * k1 - h * k2 + h * k3, u_0 + h * q1 - h * q2 + h * q3);
		q4 = f2(x + h, v_0 + h * k1 - h * k2 + h * k3, u_0 + h * q1 - h * q2 + h * q3);
		v_n += h * (k1 + 3 * k2 + 3 * k3 + k4) / 8;
		u_n += h * (q1 + 3 * q2 + 3 * q3 + q4) / 8;
	}
	return (u_n + v_n + 3.5625);
}

double find_p()
{
	double p1 = 0;
	double p2 = -1.;
	double p = (p2 - p1) / 2.;
	while (rk4(p1) > 1e-10)
	{
		p = (p2 + p1) / 2.;
		if (rk4(p1) * rk4(p) < 0)
		{
			p2 = p;
		}
		else if (rk4(p) * rk4(p2) < 0)
		{
			p1 = p;
		}
	}
	return p1;
}

void z1_diff_eq() {
	double v_next, x_next, v_prev = 0.09, x_prev = 0;
	double h;
	vector<double> n{ 10,100,1000,10000 };

	for (int i = 0; i < n.size(); i++) {
		ofstream file_to_cout_V("z1_v" + to_string(i) + ".txt");
		ofstream file_to_cout_X("z1_x" + to_string(i) + ".txt");
		v_prev = 0.09, x_prev = 0;
		v_next = 0, x_next = 0;
		h = 10. / (n[i] - 1);


		for (int j = 0; j < n[i]; j++) {
			x_next = v_prev * h + x_prev;
			v_next = U(x_next) * h + v_prev;
			file_to_cout_V << v_prev << endl;
			file_to_cout_X << x_prev << endl;
			v_prev = v_next;
			x_prev = x_next;
		}
		file_to_cout_V.close();
		file_to_cout_X.close();
	}

}

void z2_diff_eq() {

	double v_next, x_next, v_prev_prev = 0.09, x_prev_prev = 0, x_prev, v_prev;
	double h;
	vector<double> n{ 10,100,1000,10000 };

	for (int i = 0; i < n.size(); i++) {
		ofstream file_to_cout_V("z2_v" + to_string(i) + ".txt");
		ofstream file_to_cout_X("z2_x" + to_string(i) + ".txt");
		v_prev_prev = 0.09, x_prev_prev = 0;
		h = 10. / (n[i] - 1);
		x_prev = v_prev_prev * h + x_prev_prev;
		v_prev = U(x_prev) * h + v_prev_prev;

		file_to_cout_V << v_prev_prev << endl;
		file_to_cout_X << x_prev_prev << endl;
		file_to_cout_V << v_prev << endl;
		file_to_cout_X << x_prev << endl;



		for (int j = 0; j < n[i]; j++) {
			x_next = (3 * v_prev - v_prev_prev) * h / 2 + x_prev;
			v_next = (3 * U(x_next) - U(x_prev)) * h / 2 + v_prev;
			file_to_cout_V << v_next << endl;
			file_to_cout_X << x_next << endl;
			v_prev_prev = v_prev;
			v_prev = v_next;
			x_prev = x_next;
		}
		file_to_cout_V.close();
		file_to_cout_X.close();
	}



}

void z3_diff_eq() {
	double kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4;
	double h;
	double v_prev, v_next, x_prev, x_next;
	vector<double> n{ 10,100,1000,10000 };

	for (int i = 0; i < n.size(); i++) {
		ofstream file_to_cout_V("z3_v" + to_string(i) + ".txt");
		ofstream file_to_cout_X("z3_x" + to_string(i) + ".txt");
		v_prev = 0.09, x_prev = 0;
		v_next = 0, x_next = 0;
		h = 10. / (n[i] - 1);
		file_to_cout_V << v_prev << endl;
		file_to_cout_X << x_prev << endl;



		for (int j = 0; j < n[i]; j++) {
			kx1 = v_prev;
			kx2 = v_prev + h / 2. * kx1;
			kx3 = v_prev + h / 2. * kx2;
			kx4 = v_prev + h * kx3;
			kv1 = U(x_prev);
			kv2 = U(x_prev + h / 2. * kv1);
			kv3 = U(x_prev + h / 2. * kv2);
			kv4 = U(x_prev + h * kv3);
			x_next = (kx1 + 2 * kx2 + 2 * kx3 + kx4) * h / 6 + x_prev;
			v_next = (kv1 + 2 * kv2 + 2 * kv3 + kv4) * h / 6 + v_prev;
			file_to_cout_V << v_prev << endl;
			file_to_cout_X << x_prev << endl;
			v_prev = v_next;
			x_prev = x_next;
		}
		file_to_cout_V.close();
		file_to_cout_X.close();
	}

}

void z4_diff_eq() {
	double h;
	double b = 1, a = 0;
	vector<double> n{ 10,100,1000,10000 };
	for (int i = 0; i < n.size(); i++) {
		ofstream file_to_cout_V("z4_u" + to_string(i) + ".txt");
		vector<double> A(n[i]), B(n[i]), C(n[i]), F(n[i]);
		h = (b - a) / n[i];


		A[0] = 0; // edge
		B[0] = 1; // in edge indexes
		C[0] = 0; // in edge indexes
		F[0] = 0;
		for (int j = 1; j < n[i] - 1; j++) {
			A[j] = 1 / h / h + 2 / h;
			B[j] = -2 / h / h;
			C[j] = 1 / h / h - 2 / h;
			F[j] = func(a + j * h);
		}
		A[n[i] - 1] = -1 / h;
		B[n[i] - 1] = 1 + 1 / h;
		C[n[i] - 1] = 0; //edge of sys
		F[n[i] - 1] = -3.5625;

		vector<double> ANS = TripleDiag(A, B, C, F);
		for (int i = 0; i < ANS.size(); i++) {
			file_to_cout_V << ANS[i] << endl;
		}
	}
}

void z5_diff_eq() 
{
	//u''-4u'=6x^2+1
	//u(a)=0
	//u(b)+v(b)=-3.5625
	ofstream X_5("X5_out.txt");
	ofstream V_5("V5_out.txt");
	ofstream U_5("U5_out.txt");
	int n = 1000;
	double p = find_p();
	double k1, k2, k3, k4;
	double q1, q2, q3, q4;
	double h = 1. / n;
	double v_0, u_n, v_n, u_0;
	u_0 = v_0 = 0;
	u_n = p;
	v_n = u_n - 3.5625;
	double x = 0;
	for (int i = 0; i < n; i++)
	{
		X_5 << x << endl;
		U_5 << u_n << endl;
		V_5 << v_n << endl;
		k1 = f1(x, v_n, u_n);
		q1 = f2(x, v_n, u_n);
		k2 = f1(x + h / 3, v_n + h * k1 / 3, u_n + h * q1 / 3);
		q2 = f2(x + h / 3, v_n + h * k1 / 3, u_n + h * q1 / 3);
		k3 = f1(x + 2 * h / 3, v_n - h * k1 / 3, u_n - h * q1 / 3 + h * q2);
		q3 = f2(x + 2 * h / 3, v_n - h * k1 / 3, u_n - h * q1 / 3 + h * q2);
		k4 = f1(x + h, v_n + h * k1 - h * k2 + h * k3, u_n + h * q1 - h * q2 + h * q3);
		q4 = f2(x + h, v_n + h * k1 - h * k2 + h * k3, u_n + h * q1 - h * q2 + h * q3);
		v_n += h * (k1 + 3 * k2 + 3 * k3 + k4) / 8;
		u_n += h * (q1 + 3 * q2 + 3 * q3 + q4) / 8;
		x += i * h;
	}
	V_5.close();
	X_5.close();
	U_5.close();

}

int main() {
	int ZADANIE = 0;
	cout << "Ex: ";
	cin >> ZADANIE;
	switch (ZADANIE)
	{
	case 1:
		z1_diff_eq();
		break;
	case 2:
		z2_diff_eq();
		break;
	case 3:
		z3_diff_eq();
		break;
	case 4:
		z4_diff_eq();
		break;
	case 5:
		z5_diff_eq();
		break;
	}
}
