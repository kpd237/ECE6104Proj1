%=====================================
% validate_accelmove - Boris routine and particle pusher
%	Kevin Diomedi
%
%	Boris acceleration and particle movement
%=====================================
err=1e-6;


%=====================================
% Particle movement
%=====================================
%fields all zero so we have no acceleration
% just push particles with given initial velocity
theta=0;
qmr=-1;
dt=1;
lx=10;
ng=100;
E=zeros(ng,1);
B=0;
dx=lx/ng;
particle=zeros(1,4);%columns 1-xn 2-vx 3-vy 4-vz
%particle movement in middle of simulation domain
particle(1,:)=[4.5 1 5 5];%set vy, vz high to make sure they don't affect x position.
particle=accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
assert(abs(particle(1,1)-5.5)<err,'Particle Mover - Center Domain test failed');
%also test no acceleration
assert(abs(particle(1,2)-1)<err,'Particle Mover - undriven vx acceleration test failed');
assert(abs(particle(1,3)-5)<err,'Particle Mover - undriven vy acceleration test failed');
assert(abs(particle(1,4)-5)<err,'Particle Mover - undriven vz acceleration test failed');

%periodic bc on right
particle(1,:)=[lx-dx 1 5 5];%set vy, vz high to make sure they don't affect x position.
%initial speed is 1m/s, timestep is 1s.
%initial position is 9.9, final is 10.9 BC's make it 0.9.
particle=accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
assert(abs(particle(1,1)-0.9)<err,'Particle Mover - Right Periodic Boundary Condition test failed');
%periodic bc on left
particle(1,2)=-1; %set vx=-1 so we should move back to x=9.9;
particle=accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
assert(abs(particle(1,1)-9.9)<err,'Particle Mover - Left Periodic Boundary Condition test failed');


%function particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx)


%=====================================
% electric field acceleration
%=====================================
theta=0;
qmr=-1;
dt=1;
lx=10;
ng=100;
E=zeros(ng,1);
B=0;
dx=lx/ng;
particle=zeros(1,4);%columns 1-xn 2-vx 3-vy 4-vz
%test proper indexing by 
	%setting E above and below particle position
	%below
	particle(1,1)=11*dx;%intial position x
	E(11)=8;%E 1 cell below particle, should not affect it.
	particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
	assert(abs(particle(1,1)-11*dx)<err,'EField Accelerator - position 2 grid below test fail');
	assert(abs(particle(1,2)-0)<err,'EField Accelerator - vx 2 grid below test fail');
	%on - should now move particle.
	E(12)=8;
	particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
	assert(abs(particle(1,1)-3.1)<err,'EField Accelerator - grid above test fail');
	assert(abs(particle(1,2)-(-8))<err,'EField Accelerator - vx grid above test fail');

	%then test weighting by placing particle in center and setting E to be negative above and below
		%should be no acceleration
	particle(1,:)=[11.5*dx 0 0 0];
	E=zeros(ng,1);
	E(10)=-5;%cell below
	E(11)=5;%cell above
	particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
	assert(abs(particle(1,1)-dx*11.5)<err,'EField Accelerator - No movement, balanced E test fail');
	assert(abs(particle(1,2)-0)<err,'EField Accelerator - vx No movement, balanced E test fail');

	%then test with positions 25% and 75% of dx
	startx=11.25*dx;
	particle(1,:)=[startx 0 0 0];
	E=zeros(ng,1);
	E(12)=5;
	particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
	testv=qmr*dt*E(12)*0.75;
	testx=mod(testv*dt+startx,lx);
	assert(abs(particle(1,1)-testx)<err,'EField Accelerator - grid .25 test fail');
	assert(abs(particle(1,2)-testv)<err,'EField Accelerator - vx grid .25 test fail');
	%75
	E=zeros(ng,1);
	startx=11.75*dx;
	particle(1,:)=[startx 0 0 0];
	E(12)=5;
	testv=qmr*dt*E(12)*0.25;
	testx=mod(testv*dt+startx,lx);
	particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
	assert(abs(particle(1,1)-testx)<err,'EField Accelerator - grid .75 test fail');
	assert(abs(particle(1,2)-testv)<err,'EField Accelerator - vx grid .75 test fail');
	

%now test BC by performing the above at right end

%now at left end
