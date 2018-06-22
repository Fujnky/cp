//CP Blatt 7, Aufgabe 1&2. Dag-Björn Hering und Lars Funke: Euler & Runge-Kutta
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
        F << y.tail(dim * (order - 1)), f(t, y.head(dim));
        return F;
    }

    function<VectorXd(double, VectorXd)> f;
    uint order;
    uint dim;
};

//Erbende Klassen implementieren verschiedene Verfahren
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
    double alpha = 1;
    RungeKutta4ODESolver kepler([G, &alpha](double t, VectorXd r) {
      return -(alpha * G * r) / pow(r.norm(), alpha+2);
    }, 2, 3);
    MatrixXd kepler_init(2, 3);
    kepler_init << 1,   0, 0,  //     x,   y,   z
                   0, 1.1, 0;  //   v_x, v_y, v_z

    MatrixXd kepler_mars(2, 3);
    kepler_mars << -1.1,  0, 0,
                    0,  1.2, 0;

    std::ofstream kepler_file("build/rk4_kepler.txt");
    kepler_file << kepler(kepler_init, 500, 20000);

    std::ofstream kepler_file2("build/rk4_kepler2.txt");
    kepler_file2 << kepler(kepler_mars, 500, 20000);

    alpha = 1.1;
    std::ofstream kepler_alpha1("build/rk4_kepler_alpha_1.1.txt");
    kepler_alpha1 << kepler(kepler_init, 500, 20000);

    alpha = 0.9;
    std::ofstream kepler_alpha2("build/rk4_kepler_alpha_0.9.txt");
    kepler_alpha2 << kepler(kepler_init, 500, 20000);

    // d/e)
    double m1, m2, m3;
    m1 = m2 = m3 = 1;
    RungeKutta4ODESolver threebody([G,m1,m2,m3](double t, VectorXd r) {
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
    }, 2, 9);

    MatrixXd threebody_init(2, 9);
    // Die angegebenen Anfangsbedingungen haben leider nicht funktioniert.
    // Mit diesen Werten sieht es wenigstens lustig aus :)
    //                 Sonne             Planet               Mond
    threebody_init << 1, 0, 0,           0,  1, 0,            0, 0,  1.00000002,
                      1, 0, 0.0000006,  -1, -1, 0.20000001,   0, 0, -0.5;
    std::cout << threebody(threebody_init, 120, 10000);

    return EXIT_SUCCESS;
}
