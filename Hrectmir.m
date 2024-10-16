% function H defined by Eq.(29)
function H = Hrectmir(mu, nu, mup, nup, f, p)
    H = 0;
    for m=0:p.Nm
        for n=0:p.Nn
            htmp = FGrectmir(m,mu,p.Nmu,p.Dx,p.dx)*FGrectmir(m,mup,p.Nmu,p.Dx,p.dx)* ...
                   FGrectmir(n,nu,p.Nnu,p.Dy,p.dy)*FGrectmir(n,nup,p.Nnu,p.Dy,p.dy)/...
                (Kzmir(f,m,n,p)*tau(m)*tau(n));
            H = H + htmp;
        end
    end
    k    = 2*pi*f/p.c0;
    coef = k*p.Nmu*p.Dx*p.Nnu*p.Dy/((p.Dx-p.dx)*(p.Dy-p.dy));
    H    = coef*H;
end