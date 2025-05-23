clc; clear;% close all

TSDEdata1_3 = importdata("TSDEdata1_3.txt").data;
TSDEdata3 = importdata("TSDEdata3.txt").data;

x = linspace(-0.5,1.5,100)'; z = linspace(-0.5,1.5,100)';
Nx = length(x); Nz = length(z);
I = Nx; J = Nz;
M = 0.85;
beta = sqrt(1-M^2);
K = 3;
k = 0;
t = (((1-M^2)/K)^(3/2))*(M^-2);
dx = (max(x)-min(x))/Nx; dz = dx;
g = 1.4;

for i = 1:Nx
    if x(i)>0 && x(i)<1
        z(i) = 2*t*(0.25-(x(i,1)-0.5)^2);
        zprime(i) = 2*t*(1-2*x(i));
    else
        z(i) = 0;
    end
end

maxerror = 1E-6;
error = 2*maxerror;
esum = zeros(I,J);

phi = zeros(I,J);
oldphi = zeros(I,J);

while error > maxerror
    if k > 10000
        break
    end
    k = k+1;
% boundary conditions
    for j = 1:J
        phi(j,1) = phi(j,3);
    end
    for i = 2:I
        for j = 1:J
            if i == I
                phi(j,i) = phi(j,i-2);
            elseif j == 1
                if x(i) < 0 || x(i) > 1
                    phi(j,i) = phi(3,i);
                else
                    phi(j,i) = phi(3,i) - 2*dz*zprime(i);
                end
            elseif j == J
                phi(j,i) = phi(j-2,i);
            else
                phi(j,i) = (beta^2*(phi(j,i-1)+oldphi(j,i+1))/dx^2+(phi(j-1,i)+oldphi(j+1,i))/dz^2)*(2*beta^2/dx^2+2/dx^2)^-1;
            end
        end
    end
    for i = 1:I
        for j = 1:J
            esum(i,j) = abs(phi(j,i)-oldphi(j,i));
            oldphi = phi;
        end
    end
end
for i = 2:I-1
    for j = 2:J-1
        u(j,i) = (phi(j,i+1)-phi(j,i-1))/2/dx;
        cp(j,i) = -2*u(j,i)/t^(2/3);
    end
end

z_theory_3 = 2*((((1-M^2)/3)^(3/2))*(M^-2))*(0.25-(TSDEdata3(:,1)-0.5).^2);
c_d_3 = 2*trapz(z_theory_3,TSDEdata3(:,2))
z_theory_1_3 = 2*((((1-M^2)/1.3)^(3/2))*(M^-2))*(0.25-(TSDEdata1_3(:,1)-0.5).^2);
c_d_1_3 = 2*trapz(z_theory_1_3,TSDEdata1_3(:,2))
cd = 2*trapz(z(26:75),cp(2,26:75))

% figure
% surf(phi)
% figure
% surf(cp)
figure
plot(x(1:Nx-1)',cp(2,:))
hold on; axis ij; grid on
plot(TSDEdata3(:,1),TSDEdata3(:,2))
xlim([0,1])
title("Pressure Coefficient Distribution, K = 3")
xlabel("x");ylabel("$C_{p}$","interpreter","latex")
legend("Simulation","TSDE\_19","Location","northwest")

figure
plot(x,z,'k')
hold on; grid on; axis equal
plot(x,-z,'k')
xlabel("x"); ylabel("Z");
title("Airfoil Profile, K = 3")
xlim([0,1])