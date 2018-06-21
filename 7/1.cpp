//CP Blatt 7, Aufgabe 1. Dag-Björn Hering und Lars Funke: Euler & Runge-Kutta
#include <iostream>
#include <functional>
#include <Eigen/Dense>
#include <fstream>

using std::function;
using namespace Eigen;

class ODESolver {
public:
    ODESolver(const function<VectorXd(double, VectorXd)> &f, uint order, uint dim) : order(order), dim(dim) {
        shift_and_apply = [=](double t, VectorXd y) -> VectorXd {
            VectorXd F(order * dim);
            F << y.tail(dim * (order - 1)), f(t, y.head(dim));
            return F;
        };
    };

    virtual MatrixXd operator()(MatrixXd initial, double duration, uint steps) = 0;

protected:
    MatrixXd initialize(MatrixXd &initial, double duration, uint steps) {
        assert(initial.rows() == order && initial.cols() == dim);
        initial.transposeInPlace();
        initial.resize(order * dim, 1);
        MatrixXd Y = MatrixXd::Zero(dim * order, steps + 1);
        Y.col(0) = initial;
        return Y;
    }

    function<VectorXd(double, VectorXd)> shift_and_apply;
    uint order;
    uint dim;
};

class EulerODESolver : public ODESolver {
public:
    using ODESolver::ODESolver;
    MatrixXd operator()(MatrixXd initial, double duration, uint steps) override {
        auto Y = initialize(initial, duration, steps);
        double h = duration / steps;
        for (size_t i = 1; i <= steps; i++)
            Y.col(i) = h * shift_and_apply(h * i, Y.col(i - 1)) + Y.col(i - 1);
        return Y;
    };
};

class RungeKutta2ODESolver : public ODESolver {
public:
    using ODESolver::ODESolver;
    MatrixXd operator()(MatrixXd initial, double duration, uint steps) override {
        auto Y = initialize(initial, duration, steps);
        double h = duration / steps;
        for (size_t i = 1; i <= steps; i++) {
            VectorXd k_1 = h * shift_and_apply(h * (i - 1), Y.col(i - 1));
            VectorXd k_2 = h * shift_and_apply(h * (i - 1) + 0.5 * h, Y.col(i - 1) + 0.5 * k_1);
            Y.col(i) = Y.col(i - 1) + k_2;
        }
        return Y;
    };
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
    MatrixXd anfang(2, 3);
    anfang << 1, 1, 1,
              0, 0, 0;

    MatrixXd anfang_har(2, 3);
    anfang_har << 3, 2, 1,
                  0, 0, 0;

    MatrixXd anfang_np(2, 3);
    anfang_np << 1, 0, 0,
                 0, 1, 0;

    RungeKutta4ODESolver rk4([](double t, VectorXd r) { return -r; }, 2, 3);
    std::ofstream f("build/rk4.txt");
    std::ofstream o("build/rk4_har.txt");
    std::ofstream u("build/rk4_np.txt");
    f << rk4(anfang, 10, 2000) << std::endl;
    o << rk4(anfang_har, 10, 2000) << std::endl;
    u << rk4(anfang_np, 10, 2000) << std::endl;

    // Aufgabe 2
    double G = 1;
    RungeKutta4ODESolver kepler([G](double t, VectorXd r) { return -(G * r) / pow(r.norm(), 3); }, 2, 3);
    MatrixXd kepler_init(2, 3);
    kepler_init << 1, 0,   0,
                   0, 1.1, 0;

    std::ofstream kepler_file("build/rk4_kepler.txt");
    kepler_file << kepler(kepler_init, 500, 100000);
    return EXIT_SUCCESS;
}
