% function F (G) defined by Eq.(25) and Eq.(25)
function Y = FGrectmir(m, mu, N, D, d)

if (m>0)
    Y = 4*cos(m*pi*(2*mu-1)/(2*N))*sin(m*pi*(D-d)/(2*D*N))/(m*pi);
else
    Y = (D-d)/(D*N);
end