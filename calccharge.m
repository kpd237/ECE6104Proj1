function charge =calccharge(particle,qmr,dt,dx,lx,np,wp,ng)
%note this script assumes dx=1
%calculates charge from given particle species.

charge=zeros(ng,1);
q=wp^2/np/qmr/lx;
%q=wp^2/qmr*lx;%False Charge

for i=1:size(particle,1);
	xi=particle(i,1);%particle's x position
	jm1=floor(xi/dx);%multiply by dx for x_(j-1)
	%add 1 : (jm1+1) multiply by dx for x_(j)
	
	%determine index
	bel=mod(jm1,ng)+1;%must add 1 for because 1 indexing.
	abo=mod(jm1+1,ng)+1;
	

%weight below p_(j-1)=(xj-xi)
	charge(bel)=charge(bel)+((jm1+1)*dx-xi)/dx;
%weight above p_(j)=(xi-x(j-1))
	charge(abo)=charge(abo)+(xi-jm1*dx)/dx;
end
%multiply by charge of particle
charge=charge*q;

end
