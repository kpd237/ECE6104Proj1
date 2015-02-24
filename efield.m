function [E,phi]=efield(charge,dx,ng,lx,eps0,E)
spectral=true;
%==============
%FFT charge
if spectral
	phi=specpoisson(charge,dx);
else
%==========
%bc at 1st only
	bc=-2*diag(ones(ng-1,1))+diag(ones(lx-1,1),-1)+diag(ones(lx-1,1),1);
	phi=zeros(ng,1);
	phi(2:ng)=bc\charge(2:ng);
	phi=-phi/eps0*dx^2;
end
%FD electric field
E(2:lx)=-(phi(3:ng)-phi(1:lx-1))/(2*dx);
E(1)=-(phi(2)-phi(ng))/(2*dx);
E(ng)=-(phi(1)-phi(lx))/(2*dx);
end
