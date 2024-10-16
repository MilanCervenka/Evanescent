% function FG for circular elements and periodic symmetry
% defined by Eqs.(39)
function V = FGcircper(a,m,n,mu,nu, Nmu,Nnu,Deltax,Deltay)
    Lx = Nmu*Deltax;
    Ly = Nnu*Deltay;
    
    if ( (m == 0) && (n == 0) )
        V = pi*a^2/(Lx*Ly);
    else
        x0 = (mu-1/2)*Deltax;
        y0 = (nu-1/2)*Deltay;
        rho = sqrt( (m/Lx)^2 + (n/Ly)^2 );
        V = a*besselj(1, 2*pi*a*rho)/(rho*Lx*Ly);
        V = V*exp(-2i*pi*(m*x0/Lx + n*y0/Ly));
    end
end