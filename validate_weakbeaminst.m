%===============================================================
%  ECE 6104 Project 1: validating
%	Kevin Diomedi
%	3/12/2015
%
%===============================================================


theta=0;%polarization angle of magnetic field
qmre=-0.001;%charge to mass ratio electrons
qmrb=-1;%charge to mass ratio ions 
wpe=1;%electron plasma frequency
wpb=0.032;%ion plasma frequency
oce=0;%electron cyclotron freq.
ocb=0;%ion cyclotron freq
vde=0;%electron drift speed
vdb=16;%ion drift speed
L=100;%system length
ng=128;%grid points
nt=5000;%number of time steps
%nt=1;
ne=100;%number of electrons
nb=2000;%number of ions
dt=0.05;%time step
dx=L/ng;%grid spacing
Be=oce/qmre;%magnetic field efect on electrons
Bb=ocb/qmrb;%magnetic field effect on ions
k=wpe/vdb;%wave number of maximum growth rate, will be used to initial perturbation.
name='';

PIC1d3v_neutralized(ng,L,nt,dt,dx,ne,nb,vde,vdb,wpe,wpb,qmre,qmrb,k,Be,Bb,theta,name)
