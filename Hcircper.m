% function H defined by Eq.(42)
function H = Hcircper(mu, nu, mup, nup, f, p)
    H = 0;
    for m=-p.Nm:p.Nm
        for n=-p.Nn:p.Nn
            htmp = conj(FGcircper(p.a, m, n, mu , nu,  p.Nmu ,p.Nnu, p.Dx, p.Dy))*...
                        FGcircper(p.a, m, n, mup, nup, p.Nmu ,p.Nnu, p.Dx, p.Dy)/...
                   Kzper(f,m,n,p);
            H = H + htmp;
        end
    end
    k    = 2*pi*f/p.c0;
    coef = k*p.Lx*p.Ly/(pi*p.a^2);
    H    = coef*H;
end