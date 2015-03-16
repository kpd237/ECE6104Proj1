%Validate the boris mover by doing part 2 of HW#2 given in ECE 5164.
%This requires a 3D solver - so we will keep track of y,z positions in a separate vector.

v0=[0 1 0];
x0=[0 0 0];
dt=0.01;
nt=500;
qmr=1;%q=-1;m=-1;
theta=0;
B=10;
r=zeros(nt,3);%save position for plotting by time
%need to set grid parameters
ng=512;
lx=1;
dx=lx/ng;
%E field same # elements as grid.
E=zeros(1,ng);

%validation using ECE5164 HW2 2A (Fall 2014)
r= boris3dstrap(E,B,theta,qmr,dt,dx,ng,lx,x0,v0,nt);
h=figure(1);
plot3(r(:,1),r(:,2),r(:,3));
title('3D Particle (q=-1,m=-1) Orbit in Presence of Constant B_z=10 field');
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_2A.epsc2');
h=figure(2);
plot(r(:,1),r(:,2));
title('XY-Projection of Particle (q=-1,m=-1) Orbit in Presence of Constant B_z=10 field');
xlabel('X');
ylabel('Y');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_2A_proj.epsc2');

%HW2 2B - v0=[0,1,1]
r=zeros(nt,3);%save position for plotting by time
E=zeros(1,ng);
v0=[0 1 1];
x0=[0 0 0];
r= boris3dstrap(E,B,theta,qmr,dt,dx,ng,lx,x0,v0,nt);
h=figure(3);
plot3(r(:,1),r(:,2),r(:,3));
title('3D Particle (q=-1,m=-1) Orbit in Presence of Constant B_z=10 field with initial v_z');
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_3A.epsc2');
h=figure(4);
plot(r(:,1),r(:,2));
title('XY-Projection of Particle (q=-1,m=-1) Orbit in Presence of Constant B_z=10 field with initial v_z');
xlabel('X');
ylabel('Y');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_3A_proj.epsc2');
%q=+1;
r=zeros(nt,3);%save position for plotting by time
E=zeros(1,ng);
v0=[0 1 1];
x0=[0 0 0];
qmr=-1;
r= boris3dstrap(E,B,theta,qmr,dt,dx,ng,lx,x0,v0,nt);
h=figure(5);
plot3(r(:,1),r(:,2),r(:,3));
title('3D Particle (q=1,m=-1) Orbit in Presence of Constant B_z=10 field with initial v_z');
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_3A_-qmr.epsc2');
h=figure(6);
plot(r(:,1),r(:,2));
title('XY-Projection of Particle (q=1,m=-1) Orbit in Presence of Constant B_z=10 field with initial v_z');
xlabel('X');
ylabel('Y');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_3A_proj_-qmr.epsc2');

%===========================================
%hw2 -2c
%==========================================
%static electric field -Ex
E=ones(1,ng);
v0=[0 1 0];
x0=[0 0 0];
qmr=1;
r=zeros(nt,3);%save position for plotting by time
r= boris3dstrap(E,B,theta,qmr,dt,dx,ng,lx,x0,v0,nt);
h=figure(7);
plot3(r(:,1),r(:,2),r(:,3));
title('3D Particle (q=1,m=-1) Orbit in Presence of Constant B_z=10 field with E_x');
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_2c_Ex.epsc2');
h=figure(8);
plot(r(:,1),r(:,2));
title('XY-Projection of Particle (q=1,m=-1) Orbit in Presence of Constant B_z=10 field with E_x');
xlabel('X');
ylabel('Y');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_2c_proj_Ex.epsc2');
%change v0
E=ones(1,ng);
v0=[0 1 1];
x0=[0 0 0];
qmr=1;
r=zeros(nt,3);%save position for plotting by time
r= boris3dstrap(E,B,theta,qmr,dt,dx,ng,lx,x0,v0,nt);
h=figure(9);
plot3(r(:,1),r(:,2),r(:,3));
title('3D Particle (q=1,m=-1) Orbit in Presence of Constant B_z=10 field with E_x and v_z');
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_2c_ExVz.epsc2');
h=figure(10);
plot(r(:,1),r(:,2));
title('XY-Projection of Particle (q=1,m=-1) Orbit in Presence of Constant B_z=10 field with E_x and v_z');
xlabel('X');
ylabel('Y');
set(gcf, 'PaperUnit', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.75]);
saveas(h,'val_boris_2_2c_proj_ExVz.epsc2');

