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
    const Parameters& p;
    const std::map<std::string, dg::x::DVec>& ns;
    dg::x::DVec xc, yc;
    dg::x::DVec weights, tmp;
};
struct Record_static
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<void(dg::x::HVec&, Variables&, const dg::CartesianGrid2d&)> function;
    // function that generates the data points for the variable
};

// time - independent output (only called once)
std::vector<Record_static> diagnostics2d_static_list = {
    { "xc", "x-coordinate in Cartesian coordinate system",
        []( dg::x::HVec& result, Variables& v, const dg::CartesianGrid2d& grid ) {
            result = dg::evaluate( dg::cooX2d, grid);
        }
    },
    { "yc", "y-coordinate in Cartesian coordinate system",
        []( dg::x::HVec& result, Variables& v, const dg::CartesianGrid2d& grid ) {
            result = dg::evaluate( dg::cooY2d, grid);
        }
    },
    { "weights", "Gaussian Integration weights",
        []( dg::x::HVec& result, Variables& v, const dg::CartesianGrid2d& grid ) {
            result = dg::create::weights( grid);
        }
    }
};

struct Record
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<void(dg::x::DVec&, Variables&, std::string)> function;
    // function that generates the data points for the variable
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
    {"U", " 0.5 a_s mu_s N u_E^2",
        []( dg::x::DVec& result, Variables& v, std::string s) {
            dg::blas1::pointwiseDot( 0.5*v.p.a.at(s)*v.p.mu.at(s), v.ns.at(s),
                     v.rhs.uE2(), 0., result);
        }
    }
    // Maybe write out diffusion as well (should be really small though...)
};
struct Record1d
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<double( Variables&, std::string)> function;
    // function that generates the data points for the variable
};

std::vector<Record1d> diagnostics1d_s_list = {
    {"X", "Center of mass X position",
        []( Variables& v, std::string s ) {
            dg::blas1::axpby( 1., v.ns.at(s), -1., 1., v.tmp);
            double mass = dg::blas1::dot( v.weights, v.tmp);
            dg::blas1::pointwiseDot( v.tmp, v.xc, v.tmp);
            double xpos = dg::blas1::dot( v.weights, v.tmp);
            return xpos / mass;
        }
    },
    {"Y", "Center of mass Y position",
        []( Variables& v, std::string s ) {
            dg::blas1::axpby( 1., v.ns.at(s), -1., 1., v.tmp);
            double mass = dg::blas1::dot( v.weights, v.tmp);
            dg::blas1::pointwiseDot( v.tmp, v.yc, v.tmp);
            double ypos = dg::blas1::dot( v.weights, v.tmp);
            return ypos / mass;
        }
    },
    {"M1d", " N_s integral",
        []( Variables& v, std::string s) {
            return dg::blas1::dot(v.ns.at(s), v.weights);
        }
    },

    {"S1d", " a_s tau_s N ln N integral",
        []( Variables& v, std::string s) {
            double astaus = v.p.a.at(s)*v.p.tau.at(s);
            dg::blas1::transform( v.ns.at(s), v.tmp, dg::LN<double>());
            dg::blas1::pointwiseDot( astaus, v.ns.at(s), v.tmp, 0., v.tmp);
            return dg::blas1::dot( v.tmp, v.weights);

        }
    },
    {"U1d", " 0.5 a_s mu_s N u_E^2 integral",
        []( Variables& v, std::string s) {
            dg::blas1::pointwiseDot( 0.5*v.p.a.at(s)*v.p.mu.at(s), v.ns.at(s),
                     v.rhs.uE2(), 0., v.tmp);
            return dg::blas1::dot( v.tmp, v.weights);
        }
    }
};


} //namespace impurities

