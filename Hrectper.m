% function H defined by Eq.(15)
function H = Hrectper(mu, nu, mup, nup, f, p)
    H = 0;
    for m=-p.Nm:p.Nm
        for n=-p.Nn:p.Nn
            htmp = conj(FGrectper(m,mu,p.Nmu,p.Dx,p.dx))*FGrectper(m,mup,p.Nmu,p.Dx,p.dx)* ...
                   conj(FGrectper(n,nu,p.Nnu,p.Dy,p.dy))*FGrectper(n,nup,p.Nnu,p.Dy,p.dy)/...
                Kzper(f,m,n,p);
            H = H + htmp;
        end
    end
    k    = 2*pi*f/p.c0;
    coef = k*p.Nmu*p.Dx*p.Nnu*p.Dy/((p.Dx-p.dx)*(p.Dy-p.dy));
    H    = coef*H;
end