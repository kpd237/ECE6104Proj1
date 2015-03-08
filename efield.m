function [E,phi]=efield(charge,dx,ng,lx,eps0,E)
%FFT charge
phi=specpoisson(charge,dx);
%FD electric field

E(2:ng-1)=-(phi(3:ng)-phi(1:ng-2));
E(1)=-(phi(2)-phi(ng));
E(ng)=-(phi(1)-phi(ng-1));
%E(1)=0;
%E(2)=0;
E=E/(2*dx);

end
