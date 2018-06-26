//CP Blatt 7, Aufgabe 1&2. Dag-Björn Hering und Lars Funke: Doppelpendel
#include <iostream>
#include <functional>
#include <Eigen/Dense>
#include <fstream>

using std::function;
using namespace Eigen;

//Abstrakte Klasse für numerische DGL-Löser
class ODESolver {
public:
    ODESolver(const function<VectorXd(double, VectorXd)> &f, uint order, uint dim) : f(f), order(order), dim(dim) {};
    virtual MatrixXd operator()(MatrixXd initial, double duration, uint steps) = 0; //solution

protected:
    //Matrix initialisieren
    MatrixXd initialize(MatrixXd &initial, double duration, uint steps) {
        assert(initial.rows() == order && initial.cols() == dim);
        initial.transposeInPlace();
        initial.resize(order * dim, 1);
        MatrixXd Y = MatrixXd::Zero(dim * order, steps + 1);
        Y.col(0) = initial;
        return Y;
    }

    //y-Vektor verschieben und f anwenden. (Trafo auf 1. Ordnung)
    VectorXd shift_and_apply(double t, VectorXd y) {
        VectorXd F(order * dim);
        F << y.tail(dim * (order - 1)), f(t, y); //changed
        return F;
    }

    function<VectorXd(double, VectorXd)> f;
    uint order;
    uint dim;
};

class RungeKutta4ODESolver : public ODESolver {
public:
    using ODESolver::ODESolver;
    MatrixXd operator()(MatrixXd initial, double duration, uint steps) override {
        auto Y = initialize(initial, duration, steps);
        double h = duration / steps;
        for (size_t i = 1; i <= steps; i++) {
            VectorXd k_1 = h * shift_and_apply(h * i, Y.col(i - 1));
            VectorXd k_2 = h * shift_and_apply(h * i + 0.5 * h, Y.col(i - 1) + 0.5 * k_1);
            VectorXd k_3 = h * shift_and_apply(h * i + 0.5 * h, Y.col(i - 1) + 0.5 * k_2);
            VectorXd k_4 = h * shift_and_apply(h * i + h, Y.col(i - 1) + k_3);
            Y.col(i) = Y.col(i - 1) + (1 / 6.0) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
        }
        return Y;
    };
};

int main() {

    /*RungeKutta4ODESolver threebody([G, m1, m2, m3](double t, VectorXd r) {
      //Implementation eines Dreikörperproblems als 9D-System
      VectorXd r1 = r.head(3);
      VectorXd r2 = r.segment(3, 3);
      VectorXd r3 = r.tail(3);
      VectorXd result(9);
      result << -G * m2 * (r1 - r2) / pow((r1 - r2).norm(), 3)
                -G * m3 * (r1 - r3) / pow((r1 - r3).norm(), 3),
                -G * m3 * (r2 - r3) / pow((r2 - r3).norm(), 3)
                -G * m1 * (r2 - r1) / pow((r2 - r1).norm(), 3),
                -G * m1 * (r3 - r1) / pow((r3 - r1).norm(), 3)
                -G * m2 * (r3 - r2) / pow((r3 - r2).norm(), 3);
      return result;
    }, 2, 9);*/

    double m1, m2;
    m1 = m2 = 1;
    double L1, L2;
    L1 = L2 = 1;
    double g = 9.81;

    RungeKutta4ODESolver Deppelpondel([g, m1, m2, L1, L2] (double t, VectorXd theta ) {
      //Implementation eines Doppelpendel als 4D-System
      double theta1 = theta(0);
      double theta2 = theta(1);
      double theta1punkt = theta(2);
      double theta2punkt = theta(3);
      VectorXd result(2);
      double mu = m2 / (m1 + m2);
      double lambda = L1 / L2;
      double g1 = g/L1;
      double g2 = g/L2;
      double coeff = 1 / (1 - mu * pow(cos(theta2 - theta1), 2));
      result << coeff * (mu * g1 * sin(theta2) * cos(theta2-theta1) + mu * pow(theta1punkt, 2) * sin(theta2-theta1) * cos(theta2-theta1) - g1 * sin(theta1) - (mu / lambda) * pow(theta2punkt, 2) * sin(theta2 - theta1)) ,
      coeff * (g2 * sin(theta1) * cos(theta2-theta1) - mu * pow(theta2punkt, 2) * sin(theta2-theta1) * cos(theta2-theta1) - g2 * sin(theta2), - lambda * pow(theta1punkt, 2) * sin(theta2 - theta1));
      return result;
    }, 2, 2);





    MatrixXd pendulum_init_minus(2, 2);
    pendulum_init_minus << 0.1, -0.1*sqrt(2),
                     0, 0;
    MatrixXd pendulum_init_plus(2, 2);
    pendulum_init_plus << 0.1, 0.1*sqrt(2),
                          0, 0;

    double t = 10;
    uint steps_per_second = 500;
    std::cout << Deppelpondel(pendulum_init_plus, t, t*steps_per_second);

    return EXIT_SUCCESS;
}
