function charge =calccharge(particle,qmr,dt,dx,lx,np,wp,ng)
%note this script assumes dx=1
%calculates charge from given particle species.

charge=zeros(ng,1);

for i=1:size(particle,1);
	%Find grid above and below particle
	belind=floor(particle(i,1))+1;%bel(ow)ind(ex)
	abind=mod(belind,ng)+1;	 %ab(ove)ind(ex)
	%since we have 1 indexing, belind is actually the x value of the grid point above it.
	charge(belind)=charge(belind)+(belind-particle(i,1));
	charge(abind) =charge(abind) +(particle(i,1)-(belind-1));

end
%multiply by charge of particle
charge=charge*wp^2/np/qmr/dx;

end
