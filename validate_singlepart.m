theta=0;%polarization angle of magnetic field
qmre=-1;%charge to mass ratio electrons
qmri=0;%charge to mass ratio ions 
wpe=1;%electron plasma frequency
wpi=0;%ion plasma frequency
oce=-1;%electron cyclotron freq.
oci=0;%ion cyclotron freq
vde=0;%electron drift speed
vdi=0;%ion drift speed
L=4*pi/3;%system length
ng=128;%grid points
nt=3000;%number of time steps
nt=1;
ne=1;%number of electrons
ni=1;%number of ions
wuh=1.5*wpe;%upper hybrid freq
%wuh=1.4*wpe;
dt=0.2/wuh;%time step
dx=L/ng;%grid spacing
Be=oce/qmre;%magnetic field efect on electrons
Bi=oci/qmri;%magnetic field effect on ions
k=wuh/vdi%wave number of maximum growth rate, will be used to initial perturbationame='Upper Hybrid';
name='val_singlepart';
PIC1d3v(ng,L,nt,dt,dx,ne,ni,vde,vdi,wpe,wpi,qmre,qmri,k,Be,Bi,theta,name)
