function P = presr(H,d,T,z,x,t)

syms k real
w = (2*pi)/T;
rho = 1025;
g = 9.81;


P = (((H/2)*rho*g*cosh(k*(d + z))*cos((k*x) - (w*t)))/(cosh(k*d))) - (rho*g*z);
end

Fans = readtable('Froude_Krylov_Report.xlsx',"PreserveVariableNames",true); 
%assuming we put the forces in a table with first column being the time t and
%the second table being the force at time t
a = table2array(readtable('Wave_table2.xlsx',"PreserveVariableNames",true));
wvtb = a(:,1:2);
clc;
clearvars;

Dp = 320; %input("Value of depth = ") known
D = 10; %input("Value of diameter = ") known
Ln = 210; %known value of cylinder
H = 18; %input("Value of wave height = ") 
g = 9.81;
rho = 1025;
Ch = 2;
x = 0;
zser = linspace(0,-Ln,100);
tser = Fans(:,"t");
F1 = zeros(numel(tser),1);
F2 = F1;
F = F1;
f1 = zeros(100,1);
f2 = f1;
ksol = F;
%Ansys gives already gives values in discrete form so we are also taking
%values in that format
syms k real
for i = 1:length(tser)
    t = tser(i);
    
    x1 = linspace(0,5,40);
    x2 = linspace(5,10,40);
    x = linspace(0,10,40);
    for m = 1:40
        
        z = linspace(0,-210,100);
        for j = 1:numel(z)
           f1(j) = Ch*presr(H,Dp,T,z(j),x1(m),t);
           f2(j) = Ch*presr(H,Dp,T,z(j),x2(m),t);
        end
    end
   F1(i) = trapz(double(z),f1);
   F2(i) =  trapz(double(z),f2);
   F(i) = F2(i)-F1(i);
    ksol(i) = abs(vpasolve(F(i),k));
end

 k = mean(ksol);
L = (2*pi)/k;

Dl = Dp/L;
[ Dlo, ix ] = min( abs((wvtb(:,2) -Dl) )); % this gives the extrapolated value
% Dlo = wvtb(ix,2); This finds the closest value if coupled with the above
Lo = Dp/Dlo;
T1 = sqrt(1.56/Lo);
 error = T - T1;

