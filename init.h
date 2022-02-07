#pragma once

namespace impurities
{

std::map<std::string, dg::x::DVec> initial_conditions(
    const dg::x::CartesianGrid2d& grid,
    Equations<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>& eq,
    const Parameters& p,
    const dg::file::WrappedJsonValue& js)
{
    std::map<std::string, dg::x::DVec> y0;
    bool zero_pot = false;
    std::string zero_pot_species = "electrons";
    for( unsigned u=0; u<p.num_species; u++)
    {
        std::string s = p.species[u];
        dg::file::WrappedJsonValue init = js["species"][u]["init"];
        std::string type = init.get("type", "zero").asString();
        if( "const" == type)
        {
            y0[s] = dg::evaluate( dg::one, grid);
        }
        else if( "blob" == type)
        {
            double amp, sigma, posX, posY;
            amp = init["amplitude"].asDouble();
            sigma = init["sigma"].asDouble();
            posX = init["posX"].asDouble();
            posY = init["posY"].asDouble();
            double X = p.x[0] + posX*(p.x[1]-p.x[0]);
            double Y = p.y[0] + posY*(p.y[1]-p.y[0]);

            dg::Gaussian gaussian( X, Y, sigma, sigma, amp);
            dg::DVec temp = dg::evaluate( gaussian, grid);
            y0[s] = temp;
            if( p.mu.at(s) != 0 && p.tau.at(s) != 0)
            {
                std::string flr = init["flr"].asString();
                if( "none" == flr)
                    ;
                else if( "gamma" == flr)
                {
                    dg::apply( eq.gamma(s), temp, y0[s]);
                }
                else if( "gamma_inv" == flr)
                {
                    dg::apply( eq.gamma_inv(s), temp, y0[s]);
                }
            }
            dg::blas1::plus( y0[s], 1.);
        }
        else if( "wall" == type)
        {
            double wall_pos = init["posX"].asDouble();
            double wall_amp = init["amplitude"].asDouble();
            double wall_sigma = init["sigma"].asDouble();
            double X = p.x[0] + wall_pos*(p.x[1]-p.x[0]);
            dg::GaussianX wall( X, wall_sigma, wall_amp);
            dg::DVec temp = dg::evaluate( wall, grid);
            y0[s] = temp;
            if( p.mu.at(s) != 0 && p.tau.at(s) != 0)
            {
                std::string flr = init["flr"].asString();
                if( "none" == flr)
                    ;
                else if( "gamma" == flr)
                {
                    dg::apply( eq.gamma(s), temp, y0[s]);
                }
                else if( "gamma_inv" == flr)
                {
                    dg::apply( eq.gamma_inv(s), temp, y0[s]);
                }
            }
            dg::blas1::plus( y0[s], 1.);
        }
        else if( "zero_potential" == type)
        {
            if( p.mu.at(s) != 0)
                throw dg::Error( dg::Message(_ping_)<<"init type zero_potential only possible for  mu=0! You gave mu = "<<p.mu.at(s)<<"\n");
            zero_pot = true;
            zero_pot_species = s;
        }
        else
            throw dg::Error( dg::Message(_ping_)<<"init type "<<type<<" not recognized \n");


    }
    if( zero_pot)
    {
        y0[zero_pot_species] = dg::evaluate( dg::zero, grid);
        for( auto s : p.species)
        {
            if( s != zero_pot_species)
            {
                dg::x::DVec temp = y0[s];
                dg::apply( eq.gamma(s), y0[s], temp);
                dg::blas1::axpby( -p.a.at(s)/p.a.at(zero_pot_species), temp,
                        1., y0[zero_pot_species]);
            }
        }
    }
    return y0;
}

} // namespace imp
