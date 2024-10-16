% z-component of the wave vector for mirror symmetry
% defined by Eq. (21)
function kz = Kzmir(f, m, n, p)
    k  = 2*pi*f/p.c0;
    kx = m*pi/p.Lx;
    ky = n*pi/p.Ly;

    kz = sqrt(k^2 - kx^2 - ky^2);

    % the complex conjugate is calculated because of the exp(i*w*t) time convention
    kz = conj(kz);
end