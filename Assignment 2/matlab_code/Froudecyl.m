clc;
clearvars;

Dp = 320; %input("Value of depth = ");
D = 10; %input("Value of diameter = ");
H = 18; %input("Value of wave height = ");
T = 21; %input("Value of time period = ");
g = 9.81;
rho = 1025;
L = wvlength(Dp,T); 
k = (2*pi)/L;
omega = (2*pi)/T;
tser = linspace(0,T,100);
F1 = zeros(numel(tser),1);
F2 = F1;
F = F1;
Ch = 2;
%x = 0;
f1 = zeros(100,1);
f2 = f1;
%Assumption taken is Froude-krylov Force is F = rho*pi*(D^2)*accel/4
%because D/L is very small

for i = 1:length(tser)
    t = tser(i);
    
    x1 = linspace(0,5,40);
    x2 = linspace(5,10,40);
    x = linspace(0,10,40);
    for m = 1:40
        
        eta = (H/2)*cos(k*x(m) - omega*t);
        z = linspace(eta,-210,100);
        for j = 1:numel(z)
           f1(j) = Ch*pres(H,L,Dp,T,z(j),x1(m),t);
           f2(j) = Ch*pres(H,L,Dp,T,z(j),x2(m),t);
        end
    end
   F1(i) = trapz(double(z),f1);
   F2(i) =  trapz(double(z),f2);
   F(i) = F2(i)-F1(i);
end

    figure;
    hold on;
    plot(tser, F, '-','LineWidth',2);
    plot(tser, F1, '-','LineWidth',2);
    plot(tser, F2, '-','LineWidth',2);
    title('Froude Krylov force in horizontal direction');
    xlabel('t (sec)')
    ylabel('Force (N)')
    grid on;
    hold off;
 