#pragma once

#include "dg/algorithm.h"
#include "parameters.h"

namespace impurities
{

template< class Geometry, class Matrix, class Container >
struct Equations
{
    using value_type = dg::get_value_type<Container>;

    Equations( const Geometry& g, Parameters p);

    const Container& potential(std::string name ) const { return phi.at(name);}

    Helmholtz<Geometry, Matrix, Container>& gamma_inv(std::string name) {
        return m_multi_helm.at(name)[0];
    }
    const std::function<void( const Container&, Container&)>& gamma(std::string
            name) {return m_inv_helm.at(name);
    }
    const Container& uE2() const {return m_UE2;}
    const Container& real_n( std::string name) const{
        return m_real_n.at(name);
    }


    void operator()(double t, const std::map<std::string, Container>& y, std::map<std::string, Container>& yp);

private:
    void compute_phi( const std::map<std::string, Container>& y);
    void compute_psi( const Container& potential);

    Container m_binv; //magnetic field

    Container m_chi, m_omega;
    Container m_tilden, m_temp;
    Container m_phi;
    Container m_dyn, m_dxpsi, m_dypsi;
    // output relevant quantities
    Container m_uE2;
    std::map<std::string,Container> m_psi, m_real_n;
    std::vector<Container> m_multi_chi;

    dg::NestedGrids<Geometry, Matrix, Container> m_nested_grids;
    std::vector<Container> m_multi_weights;
    std::vector<dg::PCG<Container> > m_multi_pcg;

    std::map<std::string, std::vector<dg::Helmholtz<Geometry, Matrix,
        Container>>> m_multi_helm;
    std::vector<dg::Elliptic< Geometry, Matrix, Container >> m_multi_pol;

    std::map<std::function <void( const Container&, Container&)>> m_inv_helm;
    std::function< void( const Container&, Container&)> m_inv_pol;

    // initial guess
    dg::Extrapolation<Container> m_extra_phi;
    std::map<std::string, dg::Extrapolation> m_extra_n, m_extra_psi;

    // derivatives for each species
    std::map<std::string, dg::Elliptic< Geometry, Matrix, Container >> m_laplacianM;
    std::map<std::string, dg::Advection< Geometry, Matrix, Container>> m_adv;
    std::string, std::array<Matrix,2> m_centered;
    Parameters m_p;
};

template< class Geometry, class Matrix, class Container>
Equations< Geometry, Matrix, Container>::Equations( const Geometry& grid, Parameters p) :
    m_binv( evaluate( LinearX( p.kappa, 1.), grid)),
    m_chi( evaluate( dg::zero, grid )), m_omega(chi), m_UE2(chi),
    m_tilden(m_chi), m_temp(m_chi),
    m_phi(m_chi),
    m_dyn(m_chi), m_dxpsi(m_chi), m_dypsi(m_chi),
    m_extra_phi( 2, m_omega)
    m_p(p)
{
    m_nested_grids.construct( grid, m_p.num_stages);
    m_multi_chi = m_nested_girds.project(m_chi);
    for( unsigned u=0; u<m_p.num_stages; u++)
    {
        m_multi_weights[u] = dg::create::weights( m_nested_grids.grid(u));
        m_multi_pcg[u].construct( m_multi_chi[u], 10000);
    }
    for( auto s : m_p.species)
    {
        m_psi[s] = m_chi;
        m_real_n[s] = m_chi;
        m_extra_n[s] = m_extra_psi[s] = m_extra_phi;
        m_laplacianM[s].construct( grid, m_p.bcx[s], m_p.bcy[s], m_p.diff_dir);
        m_adv[s].construct( grid,  m_p.bcx[s], m_p.bcy[s]);
        // construct inversion for Gamma operators
        std::vector<std::function<void( const Container&, Container&)> > multi_inv_helm;
        for( unsigned u=0; u<m_p.num_stages; u++)
        {
            m_multi_helm[s].push_back( m_nested_grids.grid(u), m_p.bcx[s],
                    m_p.bcy[s], -0.5*m_p.tau[s]*m_p.mu[s], m_p.pol_dir);
            multi_inv_helm[u] = [
                    &op = m_multi_helm[s][u],
                    &pcg = m_multi_pcg[u],
                    &weights = m_multi_weights[u],
                    eps = m_p.eps_gamma,
                    u = u, s = s
                ] (const auto& y, auto& x)
                {
#ifdef MPI_VERSION
                    int rank;
                    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
                    unsigned number = 0;
                    dg::Timer t;
                    t.tic();
                    if( u == 0)
                        number = pcg.solve( op, x, y, 1., weights, eps, 1., 1);
                    else
                        number = pcg.solve( op, x, y, 1., weights, eps, 1., 10);
                    t.toc();
                    DG_RANK0 std::cout << "# Inverting Helmholtz for species "<<s<<" took "<<number << " iterations in "<<t.diff()<<"s\n";
                };
        }
        if( m_p.mu[s] == 0 || m_p.tau[s] == 0)
            m_inv_helm[s] = [](const auto& y, auto& x){ dg::blas1::copy( y, x);}
        else
            m_inv_helm[s] = [&, multi_inv = std::move(multi_inv_helm)]
                ( const auto& y, auto& x)
            {
                dg::nested_iterations( m_helm, x, y, multi_inv, nested_grids);
            }
    }
    std::string s = "potential";
    m_centered = {dg::create::dx( grid, m_p.bcx[s]), dg::create::dy(
                grid, m_p.bcy[s])};
    std::vector<std::function<void( const Container&, Container&)> > multi_inv_pol;
    for( unsigned u=0; u<m_p.num_stages; u++)
    {
        m_multi_pol.construct( m_nested_grids.grid(u), m_p.bcx[s], m_p.bcy[s],
                m_p.pol_dir);
        multi_inv_pol[u] = [
                &op = m_multi_pol[u],
                &pcg = m_multi_pcg[u],
                &weights = m_multi_weights[u],
                &eps = m_p.eps_pol[u],
                u = u
            ] (const auto& y, auto& x)
            {
#ifdef MPI_VERSION
                int rank;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif
                unsigned number = 0;
                dg::Timer t;
                t.tic();
                if( u == 0)
                    number = pcg.solve( op, x, y, 1., weights, eps, 1., 1);
                else
                    number = pcg.solve( op, x, y, 1., weights, eps, 1., 10);
                t.toc();
                DG_RANK0 std::cout << "# Inverting Polarisation equation took "<<number << " iterations in "<<t.diff()<<"s\n";
            };

    }
    m_inv_pol = [ &, multi_inv = std::move(multi_inv_pol) ]
        (const auto& y, auto& x)
    {
        dg::nested_iterations( m_multi_pol, x, y, multi_inv, m_nested_grids);
    }

}


//idx is impurity species one or two
template< class G, class M, class Container>
void Equations<G, M, Container>::compute_psi( double t, const Container& potential)
{
    m_multi_pol[0].variation(binv, potential, m_UE2); // u_E^2
    for( auto s : m_p.species)
    {
        m_extra_psi[s].extrapolate( t, m_psi[s]);
        dg::apply( m_inv_helm[s], potential, m_psi[s]);
        m_extra_psi[s].update( t, m_psi[s]);

        dg::blas1::axpby( 1., m_psi[s], -0.5*p.mu[s], m_UE2, m_psi[s]);   //psi  Gamma phi - 0.5 u_E^2
    }
}


template<class G, class Matrix, class Container>
void Equations<G, Matrix, Container>::compute_phi( double t, const std::vector<Container>& y)
{
    // Compute chi
    dg::blas1::copy( m_p.epsilon_D, m_chi);
    dg::blas1::pointwiseDot( m_binv, m_binv, m_temp);
    for( auto s : m_p.species)
        dg::blas1::pointwiseDot( m_p.a[s]*m_p.mu[s], y[s], m_temp, 1., m_chi);

    // update multi_pol object
    m_nested_grids.project( m_chi, m_multi_chi);
    for( unsigned u=0; u<m_p.num_stages; u++)
        m_multi_pol[u].set_chi( m_multi_chi[u]);

    //Compute rhs: \sum_s a_s Gamma_s N_s
    dg::blas1::copy( 0., m_omega);
    for( auto s : m_p.species)
    {
        // Solve Gamma N
        m_extra_n[s].extrapolate( t, m_real_n[s]);
        dg::apply( m_inv_helm[s], y[s], m_real_n[s]);
        m_extra_n[s].update( t, m_real_n[s]);
        dg::blas1::axpby( m_p.a[s], m_real_n[s], 1., m_omega);
    }

    // Solve for potential
    m_extra_phi.extrapolate( t, m_phi);
    dg::apply( m_inv_pol, m_omega, m_phi);
    extra_phi.update( t, m_phi);

    // Compute real_n (for output)
    // Technically, we only need it every output step but it is more convenient here
    for( auto s : m_p.species)
    {
        dg::blas1::pointwiseDot( m_p.a[s]*m_p.mu[s], y[s], m_temp, 0., m_chi);
        m_multi_pol[0].set_chi( m_chi);
        dg::apply( m_multi_pol[0], y[s], m_omega);
        // m_multi_pol is negative !
        dg::blas1::axpby( -1., m_omega, 1., m_real_n[s]);
    }

}

template< class G, class M, class Container>
void Equations< G, M, Container>::operator()(double t, const std::map<std::string, Container>& y, std::map<std::string, Container>& yp)
{
    // y[s] == n_s

    // compute m_phi
    compute_phi( t, y);
    // compute all m_psi
    compute_psi( t, m_phi);

    for( auto s : m_p.species)
    {
        // add advection term
        dg::blas2::symv( m_centered[0], m_psi[s], m_dxpsi);
        dg::blas2::symv( m_centered[1], m_psi[s], m_dypsi);

        // we need to subtract 1 because dg expects 0 boundary conditions
        if( m_p.bcx[s] == dg::DIR || m_p.bcy[s] == dg::DIR)
            dg::blas1::axpby( 1., y[s], -1., 1., m_tilden);
        else
            dg::blas1::copy( y[s], m_tilden);

        // ExB + Curv advection with updwind scheme
        dg::blas1::pointwiseDot( -1., m_binv, m_dypsi, m_v[0]);
        dg::blas1::pointwiseDot( +1., m_binv, m_dxpsi, m_v[1]);
        dg::blas1::plus( m_v[1], -m_p.tau[s]*m_p.kappa);
        m_adv[s].upwind( -1., m_v[0], m_v[1], m_tilden, 0., yp[s]);

        // Div ExB velocity
        dg::blas1::pointwiseDot( m_p.kappa, y[s], m_dypsi, 1., yp[s]);

        // Add diffusion
        dg::blas2::gemv( m_laplacianM[s], m_tilden, m_temp);
        dg::blas2::gemv( -m_p.nu_perp[s], m_laplacianM[s], m_temp, 1., yp[s]);
    }
}

}//namespace impurities

