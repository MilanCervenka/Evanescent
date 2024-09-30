function kz = Kzper(f, m, n, p)
    k  = 2*pi*f/p.c0;
    kx = 2*m*pi/p.Lx;
    ky = 2*n*pi/p.Ly;

    kz = sqrt(k^2 - kx^2 - ky^2);

    kz = conj(kz);
end