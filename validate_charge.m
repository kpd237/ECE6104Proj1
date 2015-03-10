%=====================================
% validate_charge.m - charge density calculation
%	Kevin Diomedi
%
%	Validate efield FD routine.
%=====================================
err=1e-6;

theta=0;
qmr=-1;
dt=1;
lx=10;
ng=100;
wp=1;
dx=lx/ng;
np=1;
particle=zeros(np,4);%columns 1-xn 2-vx 3-vy 4-vz
%test proper indexing by 
	%placing on grid point
	particle(1,1)=6*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expcharge=wp^2/np/qmr/lx;
	assert(abs(charge(7)-expcharge)<err,'Charge - indexing test failed');
	%.25
	particle(1,1)=6.25*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expbelow=wp^2/np/qmr/lx*.75;
	expabov=wp^2/np/qmr/lx*.25;
	assert(abs(charge(7)-expbelow)<err,'Charge - .25dx weighting below test failed');
	assert(abs(charge(8)-expabov)<err,'Charge - .25dx weighting above test failed');
	%.5
	particle(1,1)=6.5*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expbelow=wp^2/np/qmr/lx*.5;
	expabov=wp^2/np/qmr/lx*.5;
	assert(abs(charge(7)-expbelow)<err,'Charge - .5dx weighting below test failed');
	assert(abs(charge(8)-expabov)<err,'Charge - .5dx weighting above test failed');
	%.75
	particle(1,1)=6.75*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expbelow=wp^2/np/qmr/lx*.25;
	expabov=wp^2/np/qmr/lx*.75;
	assert(abs(charge(7)-expbelow)<err,'Charge - .75dx weighting below test failed');
	assert(abs(charge(8)-expabov)<err,'Charge - .75dx weighting above test failed');

%test right bc
	%placing on grid point
	particle(1,1)=ng*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expcharge=wp^2/np/qmr/lx;
	assert(abs(charge(1)-expcharge)<err,'Charge - indexing test failed');
	%.25
	particle(1,1)=(ng-1+0.25)*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expbelow=wp^2/np/qmr/lx*.75;
	expabov=wp^2/np/qmr/lx*.25;
	assert(abs(charge(ng)-expbelow)<err,'Charge - .25dx weighting below test failed');
	assert(abs(charge(1)-expabov)<err,'Charge - .25dx weighting above test failed');
	%.5
	particle(1,1)=(ng-1+0.5)*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expbelow=wp^2/np/qmr/lx*.5;
	expabov=wp^2/np/qmr/lx*.5;
	assert(abs(charge(ng)-expbelow)<err,'Charge - .5dx weighting below test failed');
	assert(abs(charge(1)-expabov)<err,'Charge - .5dx weighting above test failed');
	%.75
	particle(1,1)=(ng-1+0.75)*dx;
	charge=calccharge(particle,qmr,dt,dx,lx,np,wp,ng);
	expbelow=wp^2/np/qmr/lx*.25;
	expabov=wp^2/np/qmr/lx*.75;
	assert(abs(charge(ng)-expbelow)<err,'Charge - .75dx weighting below test failed');
	assert(abs(charge(1)-expabov)<err,'Charge - .75dx weighting above test failed');
