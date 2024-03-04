#include "interp.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <valarray>

// Creating type defs
using nljson = nlohmann::json;
using Vec = std::valarray<double>;

// Global variables
const double Rb = 0.0015;

double P = 1;
const double Q = 0.01;
const double R = 1;

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

// Creating a class for EKF with the iterpolation inside the class
// class EKF
// {
// private:
//     /* data */
// public:
// // Variables 
//     double q; // charge of the battery
//     double i; // battery current
//     int Ts; // sample time 
// // Methods
//     double EKF(double q, double i, double v, int Ts, interp::Interp1 f_ocv)
//     ~EKF
// };

// double EKF(double q, double i, double v, int Ts, interp::Interp1 f_ocv)
// {

// }

// EKF
// {
// }

double
columb_counting(double q, double i, int Ts)
{
    return q - (double(Ts) * i / 3600);
}

// f(x)
double
model_f(double q,double i, int Ts)
{
    return q - (double(Ts) * i / 3600);
}
// h(x)
double
model_h(double q, double i, interp::Interp1 f_ocv)
{
    return f_ocv(q) - i * Rb;
}
// fx(x)
double
model_fx(void )
{
    return 1;
}
// hx(x)
double
model_hx(double q, interp::Interp1 f_ocv)
{
    return f_ocv.der(q);
}
// EKF
double
EKF(double q, double i, double v, int Ts, interp::Interp1 f_ocv)
{
    // double P = P0;

    // Measurement update
    double Ht = model_hx(q, f_ocv);
    double Kt = (P * Ht) / (Ht * P * Ht + R);
    P = P - (Kt * Ht * P);
    q = q + (Kt * (v - (f_ocv(q) - i * Rb)));

    // Time update
    double Ft = model_fx();
    P = Ft * P * Ft + Q;
    return model_f(q, i, Ts);
}

int main()
{
    // Creating SOC and OCV vectors
    double q_max = 24;
    double q_init = 0.3 * q_max;

    // OCV data from datsheet
    Vec soc = linspace(0,1,11) * q_max;
    Vec ocv = {3.4397, 3.5299, 3.6063, 3.6418, 3.6832, 3.7546, 3.857, 3.926, 3.9978, 4.0718, 4.1798};

    // Function definition for f(x,y)
    interp::Interp1 f_ocv(soc, ocv, interp::NATURAL_CUBIC_SPLINE);

    // // Print the OCV and SOC
    // std::cout << "OCV: " << ocv << "\n";
    std::cout << "SOC: " << soc << "\n";
    std::cout << "Example: f_ocv(33%) = " << f_ocv(24 * 0.3) << "\n";
    std::cout << "Example: f_ocv'(33%) = " << f_ocv.der(24 * 0.3) << "\n";

    // Interpolate on a dense grid
    Vec soc_dense = linspace(0, 1, 200) * q_max;
    Vec ocv_dense = f_ocv(soc_dense);
    Vec docvdq_dense = f_ocv.der(soc_dense);

    // Get simulation data from json file
    std::ifstream f("simdata.json");
    nljson data = nljson::parse(f);

    Vec t_simdata = data["Simdata"]["Time"];
    Vec soc_simdata = data["Simdata"]["SOC"];
    Vec v_simdata = data["Simdata"]["Voltage"];
    Vec i_simdata = data["Simdata"]["Current"];

    // std::cout << "The size of the simulated data is " << t_simdata.size() << "\n";

// Columb counting
    int vec_len = t_simdata.size();
    auto Ts = 1; // sample tile

    Vec q_Clb_cnt(vec_len);
    q_Clb_cnt[0] = q_init;

    for (size_t ii = 1; ii < t_simdata.size(); ++ii)
        q_Clb_cnt[ii] = columb_counting(q_Clb_cnt[ii - 1], i_simdata[ii - 1], Ts);

// EKF
    Vec q_EKF(vec_len);
    q_EKF[0] = q_init;

    for (size_t ii = 1; ii < t_simdata.size(); ++ii)
        q_EKF[ii] = EKF(q_EKF[ii - 1], i_simdata[ii - 1], v_simdata[ii - 1], Ts, f_ocv);

    // std::cout << "EKF test " << EKF(q_EKF[1], i_simdata[1], v_simdata[1], Ts, f_ocv) << "\n";

    // Write output data to json file for plotting in pythonsimulatd
    nljson json{};
    // Data-sheet
    json["data"]["SOC"] = soc;
    json["data"]["OCV"] = ocv;
    // Data-sheet interpolated
    json["spline"]["SOC"] = soc_dense;
    json["spline"]["OCV"] = ocv_dense;
    json["spline"]["dOCV"] = docvdq_dense;
    // simulated data
    json["simdata"]["Time"] = t_simdata;
    json["simdata"]["SOC"] = soc_simdata;
    json["simdata"]["Voltage"] = v_simdata;
    json["simdata"]["Current"] = i_simdata;
    // columb counting
    json["columb_counting"]["SOC"] = q_Clb_cnt;
    // EKF
    json["EKF"]["SOC"] = q_EKF;

    std::ofstream nlofs("result.json", std::ios::out);
    nlofs << json;
    nlofs.close();

    return 0;
}


