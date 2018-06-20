//CP Blatt 7, Aufgabe 1. Dag-Björn Hering und Lars Funke: Euler & Runge-Kutta
#include <iostream>
#include <functional>

using std::array;
using std::function;

template<uint dim>
class ODESolver {
public:
    virtual array<double, dim> operator()(array<double, dim>) = 0;
};

template<uint dim>
class EulerODESolver : public ODESolver<dim> {
public:
    explicit EulerODESolver(function<array<double, dim>(double, array<double, dim>)> f) {
        //solve here or in seperate method?
    };

    array<double, dim> operator()(array<double, dim>) override {
        //return solved function values
    };
};

template<uint dim>
class RungeKuttaODESolver : public ODESolver<dim> {
public:
    explicit RungeKuttaODESolver(function<array<double, dim>(double, array<double, dim>)> f) {
        //solve here or in seperate method?
    };

    array<double, dim> operator()(array<double, dim>) override {
        //return solved function values
    };
};

int main() {
    function<array<double, 3>(double, array<double, 3>)> f = [](double t, array<double, 3> y){return array<double, 3>({t,t,t});};
    EulerODESolver<3> nummer1(f);
    auto result = nummer1(333);
	return EXIT_SUCCESS;
}