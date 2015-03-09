%=====================================
% validate_efield - FD Electric field solver
%	Kevin Diomedi
%
%	Validate efield FD routine.
%=====================================
ng=512;
L=2;
dx=L/ng;
x=dx*(0:ng-1);
eps0=1;
E=zeros(ng,1);

% case 1, phi=sin(pi);
charge=sin(pi*x)*pi^2;%This charge gives desired phi with spectral solver.
%this will give potential = sin(x);
[E,phi]=efield(charge,dx,ng,L,eps0,E);
h=figure(1);
plot(x,phi,x,charge/pi^2);
legend('\phi','charge');
%expected E
Eexp=-pi*cos(pi*x);
h=figure(2);
plot(x,E,x,Eexp);
legend('FFT','Expected');
title('Validation for \phi=sin(\pi*x)');
xlabel('x (m)');
ylabel('V/m');
saveas(h,'validate_FD_1.epsc2');

% case 2
N=512;
x=linspace(2*pi/N,2*pi,N);
dx=x(2)-x(1);
charge=-(sin(x).^2-cos(x)).*exp(cos(x));
phi_ana=exp(cos(x));
[E,phi]=efield(charge,dx,ng,L,eps0,E);
C=phi(1)-phi_ana(1);
figure(1);
Eexp=sin(x).*exp(cos(x));
plot(x,phi,x,phi_ana+C);
legend('\phi','expected');
figure(2);
plot(x,E,x,Eexp);
legend('FFT','Expected');
title('Validation for \phi=exp(cos(x))');
xlabel('x (m)');
ylabel('V/m');
saveas(h,'validate_FD_2.epsc2');

