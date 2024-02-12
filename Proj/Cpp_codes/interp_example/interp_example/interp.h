#pragma once
#include <iostream>
#include <memory>
#include <valarray>
#include <vector>

namespace interp {
enum InterpolationType { LINEAR, NATURAL_CUBIC_SPLINE, NOT_A_KNOT_CUBIC_SPLINE };
enum ExtrapolationType { NONE, ENDPOINTS };

using Vec = std::valarray<double>;

class Interp1 {
public:
    Interp1(){};
    Interp1(const Vec &x, const Vec &y, InterpolationType interpolation_type = NOT_A_KNOT_CUBIC_SPLINE,
            bool equidistant = true, ExtrapolationType extrap = NONE);
    Interp1(const Interp1 &p)
        : _x(p._x), _y(p._y), _x0(p._x0), _xmax(p._xmax), _delta(p._delta), n(p.n), eq(p.eq), _itype(p._itype),
          _xtype(p._xtype)
    {
        _coeffs = std::make_unique<double[]>(4 * n);
        for (size_t i = 0; i < 4 * n; i++)
            _coeffs[i] = p._coeffs[i];
    };
    ~Interp1() = default;

    Interp1 &operator=(const Interp1 &other)
    {
        _x = other._x;
        _y = other._y;
        _x0 = other._x0;
        _xmax = other._xmax;
        _delta = other._delta;
        n = other.n;
        eq = other.eq;
        _itype = other._itype;
        _xtype = other._xtype;
        _coeffs = std::make_unique<double[]>(4 * n);
        // _coeffs = new double(4 * n);
        for (size_t i = 0; i < 4 * n; i++)
            _coeffs[i] = other._coeffs[i];
        return *this;
    };

    Vec operator()(Vec &&x) const;
    Vec operator()(const Vec &x) const;
    double operator()(const double x) const;

    Vec der(const Vec &x) const;
    double der(const double x) const;
    Vec dder(const Vec &x) const;
    double dder(const double x) const;
    void print_coeffs() const;

private:
    Vec _x;
    Vec _y;
    double _x0;
    double _xmax;
    double _delta;
    size_t n;
    bool eq{true};
    InterpolationType _itype{NOT_A_KNOT_CUBIC_SPLINE};
    ExtrapolationType _xtype{NONE};
    std::unique_ptr<double[]> _coeffs;

    void compute_spline_coefficients(InterpolationType interpolation_type);
};
} // namespace interp
