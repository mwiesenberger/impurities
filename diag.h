#pragma once
// diag.h
#pragma once
#include "impurities.h"
#include "parameters.h"

namespace impurities
{

struct Variables
{
    Equations<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>& rhs;
    const dg::x::CartesianGrid2d& grid;
    const Parameters& p;
    const std::map<std::string, dg::x::DVec>& ns;
};

struct Record
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<void(dg::x::DVec&, Variables&, std::string)> function;
    // function that generates the data points for the variable
};

// time - independent output (only called once)
std::vector<Record> diagnostics2d_static_list = {
    { "xc", "x-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v, std::string s ) {
            result = dg::evaluate( dg::cooX2d, v.grid);
        }
    },
    { "yc", "y-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v, std::string s ) {
            result = dg::evaluate( dg::cooY2d, v.grid);
        }
    },
    { "weights", "Gaussian Integration weights",
        []( dg::x::DVec& result, Variables& v, std::string s ) {
            result = dg::create::weights( v.grid);
        }
    }
};

// time - dependent output (called periodically)
std::vector<Record> diagnostics2d_s_list = {
    {"n", "Real density in 2d",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            v.rhs.compute_real_n( s, v.ns.at(s), result);
        }
    },
    {"gy", "Gyro-center density in 2d",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            dg::blas1::copy(v.ns.at(s), result);
        }
    },
    {"psi", "Gyro-center potential",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            dg::blas1::copy(v.rhs.psi(s), result);
        }
    },
    {"vor", "Gyro-center vorticity",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            v.rhs.compute_lapPsi( s, result);
        }
    },
    {"S", " a_s tau_s N ln N",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            double astaus = v.p.a.at(s)*v.p.tau.at(s);
            dg::blas1::transform( v.ns.at(s), result, dg::LN<double>());
            dg::blas1::pointwiseDot( astaus, v.ns.at(s), result, 0., result);
        }
    },
    {"U", " 0.5 a_s tau_s N u_E^2",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            dg::blas1::pointwiseDot( 0.5*v.p.a.at(s)*v.p.mu.at(s), v.ns.at(s),
                     v.rhs.uE2(), 0., result);
        }
    }
    // Maybe write out diffusion as well (should be really small though...)
};


} //namespace impurities

