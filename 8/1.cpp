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
        VectorXd F(order *dim);
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
    double m1, m2;
    m1 = 1;
    m2 = 1;
    double L1, L2;
    L1 = 0.4;
    L2 = 1;
    double g = 9.81;

    RungeKutta4ODESolver Deppelpondel([g, m1, m2, L1, L2](double t, VectorXd theta) {
        //Implementation eines Doppelpendels als 4D-System
        const double theta1 = theta(0);
        const double theta2 = theta(1);
        const double theta1punkt = theta(2);
        const double theta2punkt = theta(3);
        VectorXd result(2);
        const double mu = m2 / (m1 + m2);
        const double lambda = L1 / L2;
        const double g1 = g / L1;
        const double g2 = g / L2;
        const double coeff = 1 / (1 - mu * pow(cos(theta2 - theta1), 2));
        result << coeff * (mu * g1 * sin(theta2) * cos(theta2 - theta1) +
                           mu * pow(theta1punkt, 2) * sin(theta2 - theta1) * cos(theta2 - theta1) - g1 * sin(theta1) +
                           (mu / lambda) * pow(theta2punkt, 2) * sin(theta2 - theta1)),
                coeff * (g2 * sin(theta1) * cos(theta2 - theta1) -
                         mu * pow(theta2punkt, 2) * sin(theta2 - theta1) * cos(theta2 - theta1) - g2 * sin(theta2) -
                         lambda * pow(theta1punkt, 2) * sin(theta2 - theta1));
        return result;
    }, 2, 2);

    //paar Startbedingungen
    MatrixXd pendulum_init_minus(2, 2);
    pendulum_init_minus << 0.1, -0.1 * sqrt(2), 0, 0;
    MatrixXd pendulum_init_plus(2, 2);
    pendulum_init_plus << 0.1, 0.1 * sqrt(2), 0, 0;
    MatrixXd pendulum_init_schoen(2, 2);
    pendulum_init_schoen << M_PI+0.01, M_PI, 0, 0;
    MatrixXd pendulum_init1(2, 2);
    pendulum_init1 << 0, 0, 0, 4.472;
    MatrixXd pendulum_init2(2, 2);
    pendulum_init2 << 0, 0, 0, 11.832;
    //gestört
    MatrixXd pendulum_init1_p(2, 2);
    pendulum_init1 << 0, 0, 0, 4.472 + 0.01;
    MatrixXd pendulum_init2_p(2, 2);
    pendulum_init2 << 0, 0, 0, 11.832 + 0.01;

    double t = 40;
    uint steps_per_second = 50;

    std::ofstream eins("build/eins.txt");
    eins << Deppelpondel(pendulum_init_schoen, t, t * steps_per_second);

    std::ofstream f("build/poincaré_+.txt");
    f << Deppelpondel(pendulum_init_plus, t, t * steps_per_second);

    std::ofstream ge("build/poincaré_quasi.txt");
    ge << Deppelpondel(pendulum_init1, 1000, 1000 * steps_per_second);

    std::ofstream h("build/poincaré_chaos.txt");
    h << Deppelpondel(pendulum_init2, 1000, 1000 * steps_per_second);

    std::ofstream wurst("build/poincaré_quasi_p.txt");
    wurst << Deppelpondel(pendulum_init1_p, 1000, 1000 * steps_per_second);

    std::ofstream kaese("build/poincaré_chaos_p.txt");
    kaese << Deppelpondel(pendulum_init2_p, 1000, 1000 * steps_per_second);

    return EXIT_SUCCESS;
}
