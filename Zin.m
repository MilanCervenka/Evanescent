% specific acoustic input impedance of a QWR
% defined by Eq.(46)
function z = Zin(f,rsh,fr, pars)
  % f - frequency
  % rsh - normalized resistance of the resistive sheet
  % fr - first resonance frequency of the QWR
  
  % eps, cor - to avoid zero value of the cotangent at resonance
  % (causes instability during the optimization)
  eps = 0;
  cor = 1 - 1i*eps; 

  % normalized (dimensionless) impedance
  znorm = rsh - 1i*cot(pi*f/(2*fr)/cor);

  % in physical units
  z = pars.rho0*pars.c0*znorm;
end