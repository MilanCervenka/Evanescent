% function F (G) defined by Eq.(9) and Eq.(10)
function Y = FGrectper(m, mu, N, D, d)
if (m ~= 0)
    Y = exp(-1i*m*pi*(2*mu-1)/N)*sin(m*pi*(D-d)/(D*N))/(m*pi);
else
    Y = (D-d)/(D*N);
end