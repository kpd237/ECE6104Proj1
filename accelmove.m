function particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx)
%Note this script assumes, dx=1 and integer spacing on x.

%E - electric field
%B - Magnetic Field
%particle - 1st column is position, 2nd column is vx, 3rd vy, 4th vz
%qmr - charge to mass ratio of particle species
%dt - timestep
%dx - grid spacing


tmag=(qmr*B*dt/2);
smag=2*tmag/(1+tmag^2);
%for each particle calculate acceleration
for i=1:size(particle,1);
	xi=particle(i,1);%particle's x position
	jm1=floor(xi/dx);%multiply by dx for x_(j-1)
	%add 1 : (jm1+1) multiply by dx for x_(j)
	
	%determine index
	bel=mod(jm1,ng)+1;%must add 1 for because 1 indexing.
	abo=mod(jm1+1,ng)+1;

	%Interpolate electric field.
	%since we have 1 indexing, belind is actually the x value of the grid point above it.
	Epart=E(bel)*((jm1+1)*dx-xi)/dx+E(abo)*(xi-(jm1)*dx)/dx;

	%Boris Mover
	%First half E accel
	particle(i,2)=particle(i,2)+qmr*dt*Epart/2;
	%First half roation
	v1=particle(i,2:4);
	v1(1)=v1(1)+particle(i,3)*tmag*cos(theta);
	v1(2)=v1(2)+tmag*(particle(i,4)*sin(theta)-particle(i,2)*cos(theta));
	v1(3)=v1(3)+tmag*(-particle(i,3)*sin(theta));
	%Second half rotation
	vp=particle(i,2:4);
	vp(1)=vp(1)+v1(2)*smag*cos(theta);
	vp(2)=vp(2)+smag*(v1(3)*sin(theta)-v1(1)*cos(theta));
	vp(3)=vp(3)+smag*(-v1(2)*sin(theta));
	%final half E accel
	particle(i,2:4)=vp;
	particle(i,2)=particle(i,2)+qmr*dt*Epart/2;
	

end

%prompt suggets this in separate file included with charge calculation
%I think it fits better here, unless we need to output some diagnostics before moving particles.
%easily vectorized.
particle(:,1)=mod((particle(:,1)+particle(:,2).*dt),lx);

end
