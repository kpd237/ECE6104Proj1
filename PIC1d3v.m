function PIC1d(ng,lx,nt,dt,dx,ne,nb,v0e,v0b,wpe,wpb,qme,qmb,k)
mov(1:nt)=struct('cdata',[],'colormap',[]);
eps0=1;
%column 1 is position, column 2 is velocity
backparts=zeros(ne,2);
backmass=wpe^2/qme^2*eps0/ne;
beamparts=zeros(nb,2);
beammass=wpb^2/qmb^2*eps0/nb;
%Diagnostics vs time matrices
We=zeros(nt,1);
kebeam=zeros(nt,1);
keback=zeros(nt,1);
totalq=zeros(nt,1);
totalphi=zeros(nt,1);
totalE=zeros(nt,1);

E=zeros(ng,1);
charge=zeros(ng,1);
%load initial particle position and speed
backparts(:,1)=linspace(0,lx,ne);
beamparts(:,1)=linspace(0,lx,nb);

backparts(:,2)=v0e;
beamparts(:,2)=v0b;

%perturbation to start instability
beamparts(:,1)=mod(beamparts(:,1)+0.01*sin(k*beamparts(:,1)),lx);

%calculate charge/field due to ions.
ioncharge=ones(ng,1);
ioncharge=-ioncharge*(wpb^2/qmb/nb*20+wpe^2/qme/ne*1);
for i=1:nt
	disp(i)
	
	%calculate charge from one particle spesies
	charge =calccharge(backparts,qme,dt,dx,lx,ne,wpe,ng);
	%now do for other particle species
	charge=charge+calccharge(beamparts,qmb,dt,dx,lx,nb,wpb,ng)+ioncharge;
	%calculate electric field
	[E,phi]=efield(charge,dx,ng,lx,eps0,E);
	%save diagnostics
	g=figure(1);
	subplot(3,1,1);
	scatter(backparts(:,1),backparts(:,2));
	hold all;
	scatter(beamparts(:,1),beamparts(:,2));
	hold off;
	axis([0 100 -5 30]);
	title('Phase Space');
	legend('Background','Beam');
	xlabel('x');
	ylabel('V_x');
	subplot(3,1,2);
	plot(E);
	title('Electric field');
	xlabel('x');
	ylabel('Field Strength V/m');	
	subplot(3,1,3);
	plot(charge)
	title('Charge density');
	xlabel('x');
	ylabel('Charge Density q/m^2');	
	mov(i)=getframe(g);
	%Electric Field Energy
	We(i)=eps0*sum(E.^2)/2;
	kebeam(i)=0.5*beammass*sum(beamparts(:,2).^2);
	keback(i)=0.5*backmass*sum(backparts(:,2).^2);
	totalq(i)=sum(charge);
	totalphi(i)=sum(phi);
	totalE(i)=sum(E);

	%Diagnostic screenshot
	if (i==1 || i==100 || i==2500 || i==500)
		h=figure(2);
		subplot(3,1,1);
		scatter(backparts(:,1),backparts(:,2));
		hold all;
		scatter(beamparts(:,1),beamparts(:,2));
		hold off;
		axis([0 100 -5 30]);
		title(strcat('Phase Space at t=',num2str(i*dt)));
		legend('Background','Beam');
		xlabel('x');
		ylabel('V_x');
		subplot(3,1,2);
		plot(E);
		title('Electric field');
		xlabel('x');
		ylabel('Field Strength V/m');	
		subplot(3,1,3);
		plot(charge)
		title('Charge density');
		xlabel('x');
		ylabel('Charge Density q/m^2');	
		sv=strcat('PhaseECharget',num2str(i));
		saveas(h,sv);
		saveas(h,strcat(sv,'.epsc2'));
	end
			
	%accelerate then move both particle species
	backparts=accelmove(E,backparts,qme,dt,dx,ng,lx);
	beamparts=accelmove(E,beamparts,qmb,dt,dx,ng,lx);
end

h=figure(3);
subplot(3,1,1);
plot(dt:dt:dt*nt,We);
title('Electric Field Energy Evolution');
ylabel('Energy (J)');
xlabel('time (s)');
subplot(3,1,2);
plot(dt:dt:dt*nt,keback);
title('Background Electron Kinetic Energy Evolution');
ylabel('Energy (J)');
xlabel('time (s)');
subplot(3,1,3);
plot(dt:dt:dt*nt,kebeam);
title('Beam Electron Kinetic Energy Evolution');
ylabel('Energy (J)');
xlabel('time (s)');
saveas(h,'EnergyEvolution');
saveas(h,'EnergyEvolution.epsc2');

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
plot(dt:dt:dt*nt,kebeam);
title('Total Electric Field Evolution');
ylabel('Electric Field Strength (V/m)');
xlabel('time (s)');
saveas(h,'QuantitiesEvolution');
saveas(h,'QantitiesEvolution.epsc2');
movie2avi(mov,'ece5160_Project3a.avi','compression','None');
end
