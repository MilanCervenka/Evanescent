% calculates indices mu and nu from linear index ind
% introduced by Eq. (18)
function [mu, nu] = indmunu(ind,Nmu)
    mu = mod(ind-1, Nmu)+1;
    nu = idivide(int32(ind-1), int32(Nmu))+1;
    nu = double(nu);
end