%===============================================================
%  ECE 6104 Project 1: Upper Hybrid
%	Kevin Diomedi
%	2/22/2015
%
%===============================================================


theta=0;%polarization angle of magnetic field
memir=0.01;%electron to ion mass ratio
qmre=-1;%charge to mass ratio electrons
qmri=0.01;%charge to mass ratio ions 
%qmri =1;% as given in project description gives me/mi=100!
wpe=1;%electron plasma frequency
wpi=0.1;%ion plasma frequency
oce=-1;%electron cyclotron freq.
oci=0;%ion cyclotron freq
vde=0;%electron drift speed
vdi=1;%ion drift speed
L=4*pi/3;%system length
ng=128;%grid points
nt=1000;%number of time steps
%nt=1;
ne=512;%number of electrons
ni=512;%number of ions
wuh=1.5*wpe;%upper hybrid freq
%wuh=1.4*wpe;
dt=0.2/wuh;%time step
dx=L/ng;%grid spacing
Be=oce/qmre;%magnetic field efect on electrons
Bi=oci/qmri;%magnetic field effect on ions
k=wuh/vdi%wave number of maximum growth rate, will be used to initial perturbation.
name='Upper Hybrid';
%name='FalseCharge_Upper Hybrid';
%name='largepert_Upper Hybrid';

PIC1d3v(ng,L,nt,dt,dx,ne,ni,vde,vdi,wpe,wpi,qmre,qmri,k,Be,Bi,theta,name)

theta=3*sqrt(memir);
wlh=wpi/sqrt(2);
L=10*sqrt(2)*pi;
dx=L/ng;%grid spacing
k=2*wlh/vdi;
name='Lower Hybrid'
%name='FalseCharge_Lower Hybrid'
%name='largepert_Lower Hybrid'
PIC1d3v(ng,L,nt,dt,dx,ne,ni,vde,vdi,wpe,wpi,qmre,qmri,k,Be,Bi,theta,name)
