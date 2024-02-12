#include "interp.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <valarray>

// Creating type defs
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

int main()
{
    // Creating SOC and OCV vectors
    double q_max = 24;
    Vec soc = linspace(0,100,11) * q_max;
    Vec ocv = {3.4397, 3.5299, 3.6063, 3.6418, 3.6832, 3.7546, 3.857, 3.926, 3.9978, 4.0718, 4.1798}; 

    // Function definition for f(x,y)
    interp::Interp1 f_ocv(soc, ocv, interp::NATURAL_CUBIC_SPLINE); 

    // Print the OCV and SOC 
    std::cout << "OCV: " << ocv << "\n";
    std::cout << "SOC: " << soc << "\n";

    // Interpolate on a dense grid
    Vec soc_dense = linspace(0, 100, 200) * q_max;
    Vec ocv_dense = f_ocv(soc_dense);
    Vec docvdq_dense = f_ocv.der(soc_dense);

    // Get simulation data from json file
    std::ifstream f("simdata.json");
    nljson data = nljson::parse(f);

    // Write output data to json file for plotting in python
    nljson json{};
    json["data"]["SOC"] = soc;
    json["data"]["OCV"] = ocv;
    json["spline"]["SOC"] = soc_dense;
    json["spline"]["OCV"] = ocv_dense;
    json["spline"]["dOCV"] = docvdq_dense;
    // json["simdata"]["Time"] = tsim;

    std::ofstream nlofs("result.json", std::ios::out);
    nlofs << json;
    nlofs.close();

    return 0;
}


