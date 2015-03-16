%For testing boris rotation, strap on 3d to boris routine
% single particle motion only, returns vector r of position
% versus time.
function r = boris3dstrap(E,B,theta,qmr,dt,dx,ng,lx,x0,v0,nt);

particle=[x0(1) v0];

r(1,1:3)=x0;
for t=2:nt
	particle = accelmove(E,B,theta,particle,qmr,dt,dx,ng,lx);
	r(t,:)=r(t-1,:)+particle(2:4).*dt;
	%integrating x position above eliminates the mod() operation.
	%set particlex to calculated above.
	particle(1)=r(t,1);
	%I had to do this for validation using 5164 hw2 2a
	%because the x position goes negative for a portion, mod()
	%created a huge jump to wrap around to the positive portion of
	% simulation domain.
end
