% This script calculates the absorption coefficient
% from a metasurface taking into account evanescent
% coupling among individual elements

% which figure from the article will be generated...
figure = "4(a)";
% figure = "4(b)";
% figure = "4(c)";
% figure = "4(d)";

if ( figure == "4(a)" )
    elements = "rectangular";
    symmetry = "periodic";
    fem_lumped_data = load("./data/fem_lumped_rect_periodic.txt");
    figure_title = "Rectangular elements, periodic symmetry";
end

if ( figure == "4(b)" )
    elements = "rectangular";
    symmetry = "mirror";
    fem_lumped_data = load("./data/fem_lumped_rect_mirror.txt");
    figure_title = "Rectangular elements, mirrir symmetry";
end

if ( figure == "4(c)" )
    elements = "circular";
    symmetry = "periodic";
    fem_lumped_data = load("./data/fem_lumped_circ_periodic.txt");    
    figure_title = "Circular elements, periodic symmetry";
end

if ( figure == "4(d)" )
    elements = "circular";
    symmetry = "mirror";
    fem_lumped_data = load("./data/fem_lumped_circ_mirror.txt");  
    figure_title = "Circular elements, mirror symmetry";
end

% Number of elements in the metasurface super-cell
Nmu = 4;
Nnu = 4;

% metasurface super-cell size [m]
Lx = 0.10;
Ly = 0.10;

% separation distance between the rectangular elements [m]
dx = 0.002;
dy = 0.002;

% radius of circular elements;
a = 0.011;

% dimensions of individual elements
Dx = Lx/Nmu;
Dy = Ly/Nnu;

% number of modes taken into account
% if Nm=Nn=0, the evanescent coupling is switched off
Nm = 8;
Nn = 8;   

% speed of sound [m/s]
c0 = 343;
% air ambient density [kg/m^3]
rho0 = 1.2;

% normalized resistance of the resistive sheet
rsh = 0.3;

% first resonance frequencies of the QWRs
% numbered according to Eq.(18)
fri = linspace(200, 1700, Nmu*Nnu);

% this permutes the QWRs positions randomly
% fri = fri(randperm(length(fri))); 

% calculating the cut-on frequency
if ( symmetry == "periodic" )
    fcut = c0/max(Lx, Ly);
end
if ( symmetry == "mirror" )
    fcut = c0/max(Lx, Ly)/2;
end

% structure with parameters
pars = struct('Nmu', Nmu, 'Nnu', Nnu, ...
              'Nm',  Nm,  'Nn',  Nn, ...
              'Lx',  Lx,  'Ly',  Ly, ...
              'Dx',  Dx,  'Dy',  Dy, ...
              'dx',  dx,  'dy',  dy, ...
              'a',   a, ...
              'c0',  c0, 'rho0', rho0);

% characteristic impedance          
Z0 = rho0*c0;

% frequencies for which the absorption coefficient is calculated
fmin = 1;
fmax = fcut-1;
Nf   = 200;
fi = linspace(fmin, fmax, Nf);

% starting the siopwatch
tic

for countf = 1 : length(fi)

    N = Nmu*Nnu;
    M = zeros(N, N);

    % in the periodic symmetry Eq.(16) holds and the H function
    % only depends on the difference of indices |mu-mup|, |nu-nup|
    if ( symmetry == "periodic" )
        for ii=1:Nmu
            for jj=1:Nnu
                if ( elements == "rectangular" )
                    HH(ii,jj) = Hrectper(1, 1, ii, jj, fi(countf), pars);
                end
                if ( elements == "circular" )
                    HH(ii,jj) = Hcircper(1, 1, ii, jj, fi(countf), pars);
                end
            end
        end
        
        % setting the elements of the system matrix of Eq.(17)
        for ind = 1 : N
            for indc = ind : N
                [mu,  nu ] = indmunu(ind, Nmu);
                [mup, nup] = indmunu(indc,Nmu);  
                H = HH(abs(mup-mu)+1, abs(nup-nu)+1);
                M(ind, indc) = Z0*H;
                % the system matrix is symmetric
                M(indc, ind) = M(ind, indc);
                % metasurface elements input impedances
                % system matrix diagonal
                if (ind == indc)
                    Zmunu = Zin(fi(countf), rsh, fri(ind), pars);
                    M(ind, indc) = M(ind, indc) + Zmunu;
                end
            end
        end
    end
    
    % the mirror symmetry
    if ( symmetry == "mirror" )
        for ind = 1 : N
            for indc = ind : N
                [mu,  nu ] = indmunu(ind, Nmu);
                [mup, nup] = indmunu(indc,Nmu);  
                if ( elements == "rectangular" )
                    H = Hrectmir(mu, nu, mup, nup, fi(countf), pars);
                end
                if ( elements == "circular" )
                    H = Hcircmir(mu, nu, mup, nup, fi(countf), pars);
                end
                
                M(ind, indc) = Z0*H;
                % system matrix is symmetric
                M(indc, ind) = M(ind, indc);
                % system matrix diagonal
                if (ind == indc)
                    Zmunu = Zin(fi(countf), rsh, fri(ind), pars);
                    M(ind, indc) = M(ind, indc) + Zmunu;
                end
            end
        end
    end

    % right side of Eq.(17)
    b  = 2*ones(N, 1);
    % solution of Eq.(17)
    Vi = M \ b;

    % calculating the area-weighted velocity on the metasurface
    % Eq.(33) or Eq.(45)
    if ( elements == "rectangular" )
        AvV = (Dx-dx)*(Dy-dy)*sum(Vi)/(Nmu*Dx*Nnu*Dy);
    end
    if ( elements == "circular" )
        AvV = pi*a^2*sum(Vi)/(Nmu*Dx*Nnu*Dy);
    end

    % reflection coefficient modulus - Eq.(34)
    Rkoefi(countf) = abs(1-Z0*AvV);
    % absorption coefficient
    Alfai(countf)  = 1 - Rkoefi(countf)^2;

    fprintf("%d / %d: f = %.1f Hz, A = %.3f\n", countf, length(fi), ...
        fi(countf), Alfai(countf));
end
toc

plot(fi, Alfai, 'b');
hold on;
plot(fem_lumped_data(:,1), fem_lumped_data(:,2), 'r--');
hold off;
legend("Analytic", "FEM-lumped", "Location","southeast");
title(figure_title);
xlabel("f (Hz)");
ylabel("\alpha");
axis([0 fcut 0 1.05]);



