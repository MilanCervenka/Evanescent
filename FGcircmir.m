% function FG for circular elements and mirror symmetry
% defined by Eqs.(43)
function V = FGcircmir(a,m,n,mu,nu, Nmu,Nnu,Deltax,Deltay)
    Lx = Nmu*Deltax;
    Ly = Nnu*Deltay;
    
    if ( (m == 0) && (n == 0) )
        V = pi*a^2/(Lx*Ly);
    else
        x0 = (mu-1/2)*Deltax;
        y0 = (nu-1/2)*Deltay;
        rho = sqrt( (m/Lx)^2 + (n/Ly)^2 );
        V = 2*a*tau(m)*tau(n)*besselj(1, pi*a*rho)/(rho*Lx*Ly);
        V = V*cos(pi*m*x0/Lx)*cos(pi*n*y0/Ly);
    end
end