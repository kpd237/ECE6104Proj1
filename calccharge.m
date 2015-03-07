function charge =calccharge(particle,qmr,dt,dx,lx,np,wp,ng)
%note this script assumes dx=1
%calculates charge from given particle species.

charge=zeros(ng,1);

for i=1:size(particle,1);
	%Find grid above and below particle using zero indexing.
	% this allows us to multiply by dx to get the x value at the grid point
	% also allow the use of mod() to account for periodicity.
	i
	particle(i,1)
	bel_grid=floor(particle(i,1)/dx)%grid number below (for 1 indexing need to add 1)/
	%bel_ind will always be between 0 and ng-1
	abov_grid=bel_grid+1%index 1 above this, we enforce periodicity
		% by taking mod below in the index.
	
	%weight charge to grid, recall bel/abov_grid is 0 indexing to allow dx multi
	%	add 1 for indexing matlab matrices.
	charge(bel_grid+1)=charge(bel_grid+1)+(dx*abov_grid-particle(i,1));
	charge(mod(abov_grid,ng)+1) =charge(mod(abov_grid,ng)+1) +(particle(i,1)-dx*bel_grid);
% note here we have imposed periodicity in charge index here.
end
%multiply by charge of particle
charge=charge*wp^2/np/qmr/dx/lx;

end
