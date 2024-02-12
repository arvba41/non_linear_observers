#include "interp.h"
#include <algorithm>
#include <exception>
#include <iostream>
#include <lapacke.h>
#include <memory>

#define B(i) i
#define C(i) n + i
#define D(i) 2 * n + i

namespace interp {

Interp1::Interp1(const Vec &x, const Vec &y, InterpolationType interpolation_type, bool equidistant,
                 ExtrapolationType extrap)
    : _x(x), _y(y), eq(equidistant), _itype(interpolation_type), _xtype(extrap)
{
    n = _x.size() - 1;
    _x0 = x[0];
    _xmax = x[n];
    _delta = x[1] - x[0]; // Only used in equidistant mode

    _coeffs = std::make_unique<double[]>(4 * n);
    // _coeffs = new double(4 * n);
    compute_spline_coefficients(interpolation_type);
}

void
Interp1::compute_spline_coefficients(InterpolationType interpolation_type)
{
    // double a[n + 1];
    // double h[n];
    std::unique_ptr<double[]> h = std::make_unique<double[]>(n);

    for (size_t i = 0; i < n; i++) {
        h[i] = _x[i + 1] - _x[i];
    }

    const int m = 3 * n; // Number of rows in matrix H
    // double H[m * m];
    std::unique_ptr<double[]> H = std::make_unique<double[]>(m * m);
    std::unique_ptr<double[]> rhs = std::make_unique<double[]>(m);
    std::fill(H.get(), H.get() + m * m, 0.0);
    std::fill(rhs.get(), rhs.get() + m, 0.0);

    // for (int i = 0; i < m * m; i++) {
    //     H[i] = 0.0;
    // }
    // for (int i = 0; i < m; i++) {
    //     rhs[i] = 0.0;
    // }
    size_t row, col;
    for (size_t i = 0; i < n; i++) {
        row = i;
        col = B(i);
        H[row + m * col] = h[i];
        col = C(i);
        H[row + m * col] = h[i] * h[i];
        col = D(i);
        H[row + m * col] = h[i] * h[i] * h[i];
        rhs[row] = _y[i + 1] - _y[i];
    }

    for (size_t i = 1; i < n; i++) {
        row = n + i - 1;
        col = B(i);
        H[row + m * col] = 1;
        col = B(i - 1);
        H[row + m * col] = -1;
        col = C(i - 1);
        H[row + m * col] = -2 * h[i];
        col = D(i - 1);
        H[row + m * col] = -3 * h[i - 1] * h[i - 1];

        row = 2 * n - 1 + i - 1;
        col = C(i);
        H[row + m * col] = 2;
        col = C(i - 1);
        H[row + m * col] = -2;
        col = D(i - 1);
        H[row + m * col] = -6 * h[i - 1];
    }

    switch (interpolation_type) {
    case NATURAL_CUBIC_SPLINE:
        row = 3 * n - 2;
        col = C(0);
        H[row + m * col] = 2;
        row = 3 * n - 1;
        col = C(n - 1);
        H[row + m * col] = 2;
        col = D(n - 1);
        H[row + m * col] = 6 * h[n - 1];
        break;
    case NOT_A_KNOT_CUBIC_SPLINE:
        row = 3 * n - 2;
        col = D(0);
        H[row + m * col] = 1;
        col = D(1);
        H[row + m * col] = -1;
        row = 3 * n - 1;
        col = D(n - 2);
        H[row + m * col] = 1;
        col = D(n - 1);
        H[row + m * col] = -1;
        break;
    default:
        throw std::runtime_error("Unknown interpolation type");
    }

    // Solve equations using LAPACK
    int32_t nhrs = 1;
    int32_t LDA = m;
    std::unique_ptr<int32_t[]> ipiv = std::make_unique<int32_t[]>(m);
    int32_t LDB = m;
    int32_t info;

    info = LAPACKE_dgesv(LAPACK_COL_MAJOR, m, nhrs, H.get(), LDA, ipiv.get(), rhs.get(), LDB);
    if (info != 0) {
        throw std::runtime_error("Error computing spline coefficients.");
    }

    // Store solution
    for (size_t i = 0; i < n; i++) {
        _coeffs[4 * i + 0] = _y[i];
        _coeffs[4 * i + 1] = rhs[B(i)];
        _coeffs[4 * i + 2] = rhs[C(i)];
        _coeffs[4 * i + 3] = rhs[D(i)];
    }
}

void
Interp1::print_coeffs() const
{
    for (size_t ni = 0; ni < n; ni++) {
        std::cout << ni << " (" << _x[ni] << "): " << _coeffs[4 * ni + 0] << ", " << _coeffs[4 * ni + 1] << ", "
                  << _coeffs[4 * ni + 2] << ", " << _coeffs[4 * ni + 3] << "\n";
    }
}

Vec
Interp1::operator()(const Vec &x) const
{
    Vec ret(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        ret[i] = this->operator()(x[i]);
    }
    return ret;
}

Vec
Interp1::operator()(Vec &&x) const
{
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = this->operator()(x[i]);
    }
    return x;
}

double
Interp1::operator()(const double x) const
{
    if (_xtype == NONE && (x < _x0 || x > _xmax)) {
        throw std::runtime_error("No extrapolation!");
    } else if (x <= _x0)
        return _y[0];
    else if (x >= _xmax) {
        return _y[n];
    }

    int ni;
    // double h = _delta;
    if (eq) {
        ni = int((x - _x0) / _delta);
    } else {
        auto it = std::upper_bound(std::begin(_x), std::end(_x), x);
        ni = (int)std::distance(std::begin(_x), it) - 1;
        // h = _x[ni + 1] - _x[ni];
    }

    double dx = (x - _x[ni]);
    double ret = 0.0;
    if (_itype == LINEAR) {
        throw std::runtime_error("Linear interpolation not supported");
        // ret = _y[ni] + dx / h * (_y[ni + 1] - _y[ni]);
    } else {
        ret = _coeffs[4 * ni + 0] + _coeffs[4 * ni + 1] * dx + _coeffs[4 * ni + 2] * dx * dx +
              _coeffs[4 * ni + 3] * dx * dx * dx;
    }
    return ret;
}

/**
 * @brief Interpolate the first order derivative of function
 *
 * @param x Vector of interpolation points
 * @return Vec
 */
Vec
Interp1::der(const Vec &x) const
{
    Vec ret(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        ret[i] = this->der(x[i]);
    }
    return ret;
}

/**
 * @brief Interpolate the first order derivative of function
 *
 * @param x Interpolation point
 * @return double
 */
double
Interp1::der(const double x) const
{
    if (_itype == LINEAR) {
        throw std::runtime_error("Derivatives only available for cubic spline interpolators");
    }

    int ni;
    if (x < _x0 || x > _xmax) {
        throw std::runtime_error("No extrapolation!");
    } else if (x == _xmax) {
        ni = n - 1;
    } else if (eq) {
        ni = int((x - _x0) / _delta);
    } else {
        auto it = std::upper_bound(std::begin(_x), std::end(_x), x);
        ni = (int)std::distance(std::begin(_x), it) - 1;
    }

    double dx = (x - _x[ni]);

    return _coeffs[4 * ni + 1] + 2 * _coeffs[4 * ni + 2] * dx + 3 * _coeffs[4 * ni + 3] * dx * dx;
}

/**
 * @brief Interpolate the second order derivative of function
 *
 * @param x Vector of interpolation points
 * @return Vec
 */
Vec
Interp1::dder(const Vec &x) const
{
    Vec ret(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        ret[i] = this->dder(x[i]);
    }
    return ret;
}

/**
 * @brief Interpolate the second order derivative of function
 *
 * @param x Interpolation point
 * @return double
 */
double
Interp1::dder(const double x) const
{
    if (_itype == LINEAR) {
        throw std::runtime_error("Derivatives only available for cubic spline interpolators");
    }

    int ni;
    if (x < _x0 || x > _xmax) {
        throw std::runtime_error("No extrapolation!");
    } else if (x == _xmax) {
        ni = n - 1;
    } else if (eq) {
        ni = int((x - _x0) / _delta);
    } else {
        auto it = std::upper_bound(std::begin(_x), std::end(_x), x);
        ni = (int)std::distance(std::begin(_x), it) - 1;
    }
    double dx = (x - _x[ni]);

    return 2 * _coeffs[4 * ni + 2] + 6 * _coeffs[4 * ni + 3] * dx;
}
} // namespace interp
