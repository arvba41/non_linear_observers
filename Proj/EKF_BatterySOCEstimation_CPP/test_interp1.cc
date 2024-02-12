#include "interp.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <valarray>

using nljson = nlohmann::json;
using Vec = std::valarray<double>;

Vec
linspace(double x0, double x1, size_t N)
{
    Vec y(N);

    double dx = (x1 - x0) / (N - 1);
    y[0] = x0;
    for (size_t k = 1; k < N - 1; ++k)
        y[k] = y[k - 1] + dx;
    y[N - 1] = x1;

    return y;
}

std::ostream &
operator<<(std::ostream &os, const Vec &v)
{
    os << "(";
    for (size_t i = 0; i < v.size() - 1; i++)
        os << v[i] << ", ";
    os << v[v.size() - 1] << ")";
    return os;
}

Vec
fun(const Vec &x)
{
    // return 1. / (1. + std::exp(-x));
    return std::sin(x);
}

int
main()
{
    // Create (x, y) pairs on a sparse grid and create interpolation object
    Vec x = linspace(-5., 5., 10);
    Vec y = fun(x);

    // Function definition for f(x,y)
    interp::Interp1 f(x, y, interp::NATURAL_CUBIC_SPLINE); 

    std::cout << "x: " << x << "\n";
    std::cout << "y: " << y << "\n";

    std::cout << "Example: f(pi / 2) = " << f(M_PI / 2) << "\n";
    std::cout << "Example: f'(pi / 2) = " << f.der(M_PI / 2) << "\n";

    // Interpolate on a dense grid
    Vec x_dense = linspace(-5, 5, 200);
    Vec y_dense = f(x_dense);
    Vec dydx_dense = f.der(x_dense);

    // Write output data to json file for plotting in python
    nljson json{};
    json["data"]["x"] = x;
    json["data"]["y"] = y;
    json["spline"]["x"] = x_dense;
    json["spline"]["y"] = y_dense;
    json["spline"]["dy"] = dydx_dense;

    std::ofstream nlofs("result.json", std::ios::out);
    nlofs << json;
    nlofs.close();

    return 0;
}
