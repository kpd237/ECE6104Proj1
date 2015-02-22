function [E,phi]=efield(charge,dx,ng,lx,eps0,E)

%==============
%FFT charge
%spect_chg=fft(charge);
%calculate potential
%kn=2*pi/ng*linspace(-ng/2,ng/2,ng);
%Kn2=(kn.^2.*(sin(kn.*dx/2)./(kn.*dx/2)).^2).';
%spect_pot=spect_chg./(eps0*Kn2);

%IFFT potential
%phi=ifft(spect_pot,ng);
%phi=real(phi);

%==========
%-2 on main diagonal version
%bc=-2*diag(ones(ng-2,1))+diag(ones(lx-2,1),-1)+diag(ones(lx-2,1),1);
%phi=zeros(ng,1);
%phi(2:lx)=bc\charge(2:lx);
%phi=-phi/eps0*dx^2
%bc at 1st only
bc=-2*diag(ones(ng-1,1))+diag(ones(lx-1,1),-1)+diag(ones(lx-1,1),1);
phi=zeros(ng,1);
phi(2:ng)=bc\charge(2:ng);
phi=-phi/eps0*dx^2;

%FD electric field
E(2:lx)=-(phi(3:ng)-phi(1:lx-1))/(2*dx);
E(1)=-(phi(2)-phi(ng))/(2*dx);
E(ng)=-(phi(1)-phi(lx))/(2*dx);
end
