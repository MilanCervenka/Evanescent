% function H defined by Eq.(44)
function H = Hcircmir(mu, nu, mup, nup, f, p)
    H = 0;
    for m=0:p.Nm
        for n=0:p.Nn
            htmp = FGcircmir(p.a, m, n, mu , nu,  p.Nmu ,p.Nnu, p.Dx, p.Dy)*...
                   FGcircmir(p.a, m, n, mup, nup, p.Nmu ,p.Nnu, p.Dx, p.Dy)/...
                   (tau(m)*tau(n)*Kzmir(f,m,n,p));
            H = H + htmp;
        end
    end
    k    = 2*pi*f/p.c0;
    coef = k*p.Lx*p.Ly/(pi*p.a^2);
    H    = coef*H;
end