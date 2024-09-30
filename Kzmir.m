function kz = Kzmir(f, m, n, p)
    k  = 2*pi*f/p.c0;
    kx = m*pi/p.Lx;
    ky = n*pi/p.Ly;

    kz = sqrt(k^2 - kx^2 - ky^2);
    kz = conj(kz);
end