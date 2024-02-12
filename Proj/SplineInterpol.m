classdef SplineInterpol < handle
    properties
        pp;
        der_pp;
    end
    methods
        function obj = SplineInterpol(x, y)
            obj.pp = spline(x, y);
            obj.der_pp = der_spline(obj.pp);
        end
        function r = f(obj, xi)
            r = ppval(obj.pp, xi);
        end
        function r = df(obj, xi)
            r = ppval(obj.der_pp, xi);
        end
    end
end

function pp_der = der_spline(pp)
    pp_der = pp;
    pp_der.coefs = pp_der.coefs*diag([3, 2, 1, 0]);
    pp_der.coefs = [zeros(size(pp_der.coefs, 1), 1) pp_der.coefs(:, 1:3)];
end
