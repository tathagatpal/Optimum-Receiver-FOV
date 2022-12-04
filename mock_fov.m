clc; clear all;

fov_vec = linspace(deg2rad(15), deg2rad(45), 10);
h0 = 3;
R = 2;
semiang = deg2rad(60);
m = -log(2)/log(cos(semiang));
A = 10^-4;
n = 1.5;


x0 = 0;


%r0 = sqrt(xx.^2 + yy.^2);
%theta_vec =atan(yy./xx);

beta_vec = linspace(0,0,1);
N = 4;
P = 2;

r0 = linspace(0, eps, 1000);
theta_vec = linspace(0,2*pi,1000);

for tb = 1:length(beta_vec)
    beta = beta_vec(tb);
    %beta = 0;
for i = 1: length(fov_vec)
%    h0 = h0_vec(i);
    for j=1:length(theta_vec)
       

        theta = theta_vec(j) + (2*pi/N)*[0:1:N-1];
               
        
        fov = fov_vec(i);
         upper = h0*tan(fov);

        for q = 1:N
            rN(q,:) = sqrt(x0^2 + r0.^2 - 2*r0*x0*cos(theta(q)));
            func(q,:) = exp(-beta*rN(q,:)).*(A*(m+1)*n^2*((h0./sqrt(h0^2+rN(q,:).^2)).^(3))./((h0^2+rN(q,:).^2)*(sin(2*fov))^2)).*((rN(q,:) < upper));
            contributionN(q) = (1/R)*trapz(r0,func(q,:));
        end

        simul(j) = sum(contributionN);



% 
% % 
%         r1 = sqrt(x0^2 + r0.^2 - 2*r0*x0*cos(theta1));
%         r2 = sqrt(x0^2 + r0.^2 - 2*r0*x0*cos(theta2));
%         r3 = sqrt(x0^2 + r0.^2 - 2*r0*x0*cos(theta3));
%         r4 = sqrt(x0^2 + r0.^2 - 2*r0*x0*cos(theta4));
%         func1 = exp(-beta*r1).*(A*(m+1)*n^2*((h0./sqrt(h0^2+r1.^2)).^(3))./((h0^2+r1.^2)*(sin(2*fov))^2));
%         func1(r1 > upper) = 0;
%         func2 = exp(-beta*r2).*(A*(m+1)*n^2*((h0./sqrt(h0^2+r2.^2)).^(3))./((h0^2+r2.^2)*(sin(2*fov))^2));
%         func2(r2 > upper) = 0;
%         func3 = exp(-beta*r3).*(A*(m+1)*n^2*((h0./sqrt(h0^2+r3.^2)).^(3))./((h0^2+r3.^2)*(sin(2*fov))^2));
%         func3(r3 > upper) = 0;
%         func4 = exp(-beta*r4).*(A*(m+1)*n^2*((h0./sqrt(h0^2+r4.^2)).^(3))./((h0^2+r4.^2)*(sin(2*fov))^2));
%         func4(r4 > upper) = 0;

       %simul(j) = (1/R) * (trapz(r0,func1) + trapz(r0,func2) + trapz(r0,func3) + trapz(r0,func4));
%        simul(j) = (2/R^2) * (trapz(r0,func1.*r0) + trapz(r0,func2.*r0) + trapz(r0,func3.*r0) + trapz(r0,func4.*r0));

%         
%         ana(i,j) = (1/(3*h0^4)) * (2*(h0*tan(fov))^3 + 3*h0*tan(fov))/((h0*sec(fov))^3 * (sin(2 * fov))^2);
%         ana2(i) = (1/(3)) * (2*R^3 + 3*R)/((R^2 + 1)^1.5 * (sin(2 * fov))^2);
    end
    final(i) = (1/(2*pi)) * trapz(theta_vec, simul);
end

[val_fov(tb) ind_fov(tb)] = max(final);

rad2deg(2*fov_vec(ind_fov))

tb
plot(rad2deg(2*fov_vec),pow2db(1e3*final))
hold on
end

%plot(beta_vec, rad2deg(2*fov_vec(ind_fov)))


%plot(rad2deg(fov_vec),ana,'*')
