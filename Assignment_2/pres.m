function P = pres(H,L,d,T,z,x,t)

k = (2*pi)/L;
w = (2*pi)/T;
rho = 1025;
g = 9.81;


P = (((H/2)*rho*g*cosh(k*(d + z))*cos((k*x) - (w*t)))/(cosh(k*d))) - (rho*g*z);
end