%========================================================
%	PIC1d3v.m
%		Kevin Diomedi
%
%	
%	Performs 1D spacial 3D velocity Particle in cell 
%		Simulation on 2 particle species assuming
%		periodic boundary conditions.
%
%========================================================

function PIC1d3v(ng,lx,nt,dt,dx,n1,n2,vd1,vd2,wp1,wp2,qm1,qm2,k,B1,B2,theta,name)
%	ng = number of grid cells
%	lx = length of system
%	nt = number of time steps
%	dt = size of time step
%	dx = size of grid cell
%	n1 = number of particles in species 1
%	n2 = number of particles in species 2
%	vd1 = initial drift velocity of species 1
%	vd2 = initial drift velocity of species 2
%	wp1 = plasma frequency of species 1
%	wp2 = plasma frequency of species 2
%	qm1 = charge to mass ratio of species 1
%	qm2 = charge to mass ratio of species 2
%	k = wavenumber of expected max grow rate, used to initial perturbation.
%	B1 = magnitude of impressed static magnetic field on species 1
%	B2 = magnitude of impressed static magnetic field on species 2
%	theta = polarization angle of impressed static magnetic field
debug=true;
debug=false;
%mov(1:nt)=struct('cdata',[],'colormap',[]);
eps0=1;
%column 1 is position, column 2 is velocity
species1=zeros(n1,4);
spec1mass=wp1^2/qm1^2*eps0/n1;
species2=zeros(n2,4);
spec2mass=wp2^2/qm2^2*eps0/n2;
%Diagnostics vs time matrices
We=zeros(nt,1);
kespec2=zeros(nt,1);
kespec1=zeros(nt,1);
totalq=zeros(nt,1);
totalphi=zeros(nt,1);
totalE=zeros(nt,1);

E=zeros(ng,1);
charge=zeros(ng,1);
%load initial particle position and speed
dp1=lx/n1;%particle spacing
dp2=lx/n2;
ind=0:n1-1;
species1(:,1)=dp1*ind;
ind=0:1:n2-1;
species2(:,1)=dp2*ind;

species1(:,2)=vd1;%assign drift vx velocity to particle species1
species2(:,2)=vd2;%assign drift vx velocity to particle species2
x=(0:1:ng-1)*dx;

%perturbation to start instability
perbx=lx/n1/2;
perbv=perbx/2;
%perbv=0.0;
%perbx=0.5;%large pert
%perbv=0.25;%large pert
species1(:,2)=species1(:,2)+perbv*sin(k*species1(:,1));%Vx perturbation
%species1(:,3)=species1(:,3)+perbv*cos(k*species1(:,1));%Vy perturbation
%species1(:,4)=species1(:,4)+perbv*cos(k*species1(:,1));%Vz perturbation
species1(:,1)=mod(species1(:,1)+perbx*cos(k*species1(:,1)),lx);%x peturbation
%species2(:,1)=mod(species2(:,1)+perbx*cos(k*species2(:,1)),lx);

for i=1:nt
	disp(i)
	
	%accelerate then move both particle species
	species1=accelmove(E,B1,theta,species1,qm1,dt,dx,ng,lx);
	species2=accelmove(E,B2,theta,species2,qm2,dt,dx,ng,lx);

	%calculate charge from one particle spesies
	charge =calccharge(species1,qm1,dt,dx,lx,n1,wp1,ng)+calccharge(species2,qm2,dt,dx,lx,n2,wp2,ng);
	%calculate electric field
	[E,phi]=efield(charge,dx,ng,lx,eps0,E);
	if debug
		%For these plots to come out correctly
		% Make sure k=1.5
		% no perturbation on species 1 (electrons)
		% 0.01 perturbation on species 2 (ions)
		% comment out the accelmov above
		g=figure(15);
		subplot(2,1,1);
		plot(x,E,x,-1/1.5*0.01/512*cos(1.5*x));
		title('Electric field validation 0.01sin(1.4*x) perturbation');
		xlabel('x');
		ylabel('Field Strength V/m');	
		subplot(2,1,2);
		plot(x,charge,x,0.01/512*sin(1.5*x))
		title('Charge density validation 0.01sin(1.4*x) perturbation');
		xlabel('x');
		ylabel('Charge Density q/m^2');	
		legend('Simulation','Expected');
		saveas(g,strcat(name,'validation_initial_analytical.epsc2'));
		
	end
	%save diagnostics
	g=figure(1);
	subplot(4,1,1);
	scatter(species1(:,1),species1(:,2));
	hold all;
	scatter(species2(:,1),species2(:,2));
	hold off;
	axis([0 lx -0.5 1.5]);
	title('Phase Space');
	legend('Electrons','Ions');
	xlabel('x');
	ylabel('V_x');
	subplot(4,1,2);
	plot(x,E);
	title('Electric field');
	xlabel('x');
	ylabel('Field Strength V/m');	
	subplot(4,1,3);
	plot(x,charge)
	title('Charge density');
	xlabel('x');
	ylabel('Charge Density q/m^2');	
	subplot(4,1,4);
	scatter(species1(:,2),species1(:,3));
	hold all;
	scatter(species2(:,2),species2(:,3));
	hold off;
	title('Velocity Space');
	xlabel('V_x');
	ylabel('V_y');	
	%mov(i)=getframe(g);
	%Electric Field Energy
	We(i)=eps0*sum(E.^2)/2;
	kespec2(i)=0.5*spec2mass*sum(species2(:,2).^2);
	kespec1(i)=0.5*spec1mass*sum(species1(:,2).^2);
	totalq(i)=sum(charge);
	totalphi(i)=sum(phi);
	totalE(i)=sum(E);

	%Diagnostic screenshot
	if (i==1 || i==100 || i==nt || i==500)
		h=figure(2);
		subplot(3,1,1);
		scatter(species1(:,1),species1(:,2));
		hold all;
		scatter(species2(:,1),species2(:,2));
		hold off;
		axis([0 lx -0.5 1.5]);
		title(strcat('Phase Space at t=',num2str(i*dt)));
		legend('Electrons','Ions');
		xlabel('x');
		ylabel('V_x');
		subplot(3,1,2);
		plot(x,E);
		title('Electric field');
		xlabel('x');
		ylabel('Field Strength V/m');	
		subplot(3,1,3);
		plot(x,charge)
		title('Charge density');
		xlabel('x');
		ylabel('Charge Density q/m^2');	
		sv=strcat(name,'PhaseECharget',num2str(i));
		saveas(h,sv);
		saveas(h,strcat(sv,'.epsc2'));
	end
			
end
h=figure(10);
plot(dt:dt:dt*nt,10*log10(We));
title('Electric Field Energy Evolution');
ylabel('Energy (dB)');
xlabel('time (s)');
saveas(h,strcat(name,'Efieldgrowth'));
saveas(h,strcat(name,'Efieldgrowth.epsc2'));

h=figure(3);
subplot(3,1,1);
plot(dt:dt:dt*nt,We);
title('Electric Field Energy Evolution');
ylabel('Energy (J)');
xlabel('time (s)');
subplot(3,1,2);
plot(dt:dt:dt*nt,kespec1);
title('Background Electron Kinetic Energy Evolution');
ylabel('Energy (J)');
xlabel('time (s)');
subplot(3,1,3);
plot(dt:dt:dt*nt,kespec2);
title('Beam Electron Kinetic Energy Evolution');
ylabel('Energy (J)');
xlabel('time (s)');
saveas(h,strcat(name,'EnergyEvolution'));
saveas(h,strcat(name,'EnergyEvolution.epsc2'));

h=figure(4);
subplot(3,1,1);
plot(dt:dt:dt*nt,totalq);
title('Total Charge Density Evolution');
ylabel('Charge Density (Q/m)');
xlabel('time (s)');
subplot(3,1,2);
plot(dt:dt:dt*nt,totalphi);
title('Total Potential Evolution');
ylabel('Potential (V)');
xlabel('time (s)');
subplot(3,1,3);
plot(dt:dt:dt*nt,kespec2);
title('Total Electric Field Evolution');
ylabel('Electric Field Strength (V/m)');
xlabel('time (s)');
saveas(h,strcat(name,'QuantitiesEvolution'));
saveas(h,strcat(name,'QantitiesEvolution.epsc2'));
%movie2avi(mov,'ece5160_Project3a.avi','compression','None');
end
