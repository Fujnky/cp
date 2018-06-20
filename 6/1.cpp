//CP Blatt 6, Aufgabe 1. Dag-Björn Hering und Lars Funke: Lanczos et al.

#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <complex>
#include <string>

#include <chrono>

using namespace Eigen;
using std::cout;
using std::endl;

MatrixXd householder(MatrixXd A){
	//Siehe letzte Abgabe von Jasper und Henning ;-)
	int N = A.rows();
	MatrixXd Arueck(N, N);
	Arueck = A;

	for (size_t n = 1; n <= N-2; n++) {//Transformation N-2 mal
		MatrixXd v(N-n, 1);
		v = Arueck.col(n-1).tail(N-n);
		double k = sqrt(v.squaredNorm()); // die quadratische Summe der letzten  (N-n) Zahlen von der n-ten Spalte
		k = copysign(abs(k), -Arueck(n, n-1)); // um Rundungsfehler zu vermeiden
		MatrixXd u(N-n, 1);
		u = (v-k*VectorXd::Unit(N-n, 0))/sqrt((v-k*VectorXd::Unit(N-n, 0)).squaredNorm());

		MatrixXd S(N-n, N-n);
		S = MatrixXd::Identity(N-n, N-n) - 2*u*u.transpose();

		MatrixXd P(N, N);
		P << MatrixXd::Identity(n, n), MatrixXd::Zero(n, N-n), MatrixXd::Zero(N-n, n),  S;
		Arueck = P.transpose()*Arueck*P;
	}
	return Arueck;
}


MatrixXd Lanczos(MatrixXd A, MatrixXd q){
	//Lanczos algorithm
	uint dim = A.rows();
	MatrixXd Q(dim, dim);
	Q.col(0) = q / q.norm(); //initialize first column
  MatrixXd T = MatrixXd::Zero(dim, dim);
	double alpha = Q.col(0).adjoint() * A * Q.col(0);
	VectorXd x = A * Q.col(0) - alpha * Q.col(0);
	T(0, 0) = alpha;
	double beta = x.norm();
	Q.col(1) = x / beta;
	T(1, 0) = T(0, 1) = beta;

	for (size_t i = 2; i < dim; i++) {
		  T(i-1, i-1) = Q.col(i-1).adjoint() * A * Q.col(i-1); //alpha
		 	VectorXd x = A * Q.col(i-1) - T(i-2, i-1) * Q.col(i-2) - T(i-1, i-1) * Q.col(i-1);
			T(i-1, i) = T(i, i-1) = x.norm(); //beta
		 	Q.col(i) = x / T(i-1, i); // next q_i
  }

	T(dim-1, dim-1) = Q.col(dim - 1).adjoint() * A *Q.col(dim - 1);
	return T;
}

MatrixXd QR_decompose(MatrixXd A){
  MatrixXd Q = A; //just for same dimension
  MatrixXd U = A;
  for (size_t k = 0; k < U.rows(); k++){
    for(size_t j = 0; j < k; j++)
			//generate u_k
		  U.col(k) -= (double)(U.col(j).adjoint() * A.col(k)) / U.col(j).squaredNorm() * U.col(j);
    Q.col(k) = U.col(k) / U.col(k).norm();
  }
  return Q;
}

MatrixXd QR_iteration(MatrixXd A, size_t steps)
{
	//iterate QR decomposition steps
	for(size_t i = 0 ; i < steps ; i++)
	{
  	auto Q = QR_decompose(A);
   	A = Q.adjoint() * A * Q;
	}
	return A;
}

double time_diff(std::chrono::high_resolution_clock::time_point start)
{
	return  std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
}

int main() {

	MatrixXd M(6, 6);
	M <<  9, 10,  -7, 15,  13, 14,
			 10,  1,   7, 15,   6, 11,
			 -7,  7,  17,  3, -16, -4,
			 15, 15,   3, -4,   7,  9,
			 13,  6, -16,  7,  11, 16,
			 14, 11,  -4,  9,  16,  8;

	MatrixXd save_time_lanczos(18,50);
	MatrixXd save_time_house(18,50);
	MatrixXd save_time_eigen(18,50);

	//time different sized matrices
	for(size_t N = 3 ; N<=20 ; N++)
		for(size_t i = 0 ; i< 50 ; i++){
			MatrixXd random = MatrixXd::Random(N,N);
			VectorXd q = VectorXd::Constant(N, 1);
			// time for Lanczos
			auto start_lanczos = std::chrono::high_resolution_clock::now();
			auto eigen_lanczos = QR_iteration(Lanczos(random, q),30).diagonal();
			save_time_lanczos(N-3,i) = time_diff(start_lanczos);
			// time for Householder
			auto start_householder = std::chrono::high_resolution_clock::now();
			auto eigen_householder = QR_iteration(householder(random),30).diagonal();
			save_time_house(N-3,i) = time_diff(start_householder);
			// time for Eigen
			auto start_eigen = std::chrono::high_resolution_clock::now();
			EigenSolver<MatrixXd> esM(random);
			auto eigen_solver = esM.eigenvalues();
			save_time_eigen(N-3,i) = time_diff(start_eigen);
	}

	std::ofstream f1("build/lanczos.txt");
	std::ofstream f2("build/householder.txt");
	std::ofstream f3("build/eigen.txt");
	f1 << save_time_lanczos.rowwise().mean();
	f2 << save_time_house.rowwise().mean();
	f3 << save_time_eigen.rowwise().mean();

	VectorXd q_0_1(6, 1);
	q_0_1 << 1, 1, 1, 1, 1, 1;

	VectorXd q_0_2(6, 1);
	q_0_2 << 1, 0, 0, 0, 0, 0;

	auto L1 = Lanczos(M, q_0_1);
	auto L2 = Lanczos(M, q_0_2);
	cout << "Lanczos 111111"<< endl << L1 << endl;
	cout << "Lanczos 100000"<< endl << L2 << endl;
	cout << "eigenvalues lanczos111111+qr:" << endl << QR_iteration(L1, 40).diagonal() << endl;
	cout << "eigenvalues lanczos100000+qr:" << endl << QR_iteration(L2, 40).diagonal() << endl;

  EigenSolver<MatrixXd> esM(M);
	cout << "eigenvalues eigen:" << endl << esM.eigenvalues() << endl;

	MatrixXd L = Lanczos(M, q_0_1);
	EigenSolver<MatrixXd> esL(L);
	cout << "Eigenwerte von L sind: " << endl << esL.eigenvalues() << endl;
	cout << "householder" << endl<< householder(M) << endl;
	return EXIT_SUCCESS;
}
