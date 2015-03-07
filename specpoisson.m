function phi = specpoisson(rho,dx)

eps0=1;

N=length(rho);
ng2=N/2;
for k=1:ng2
        kdx2=(pi/N)*k;
        ksqi(k)=+eps0/((2.0*sin(kdx2)/dx).^2);
end

rho(1)=0;
%rhok(N+1)=0;%not sure on this one)
%hdx=0.5*dx;

%rho=rho*hdx;
%Take fft and reorder.
rhok=fft(rho,pow2(nextpow2(length(rho))));
rhok=fftshift(rhok);
rhok(1)=0;

phik(1)=0;
for k=2:ng2
	kk=N+2-k;
	phik(k)=ksqi(k-1)*rhok(k);
	phik(kk)=ksqi(k-1)*rhok(kk);
	phik(ng2+1)=ksqi(ng2)*rhok(ng2+1);
end

%li=1.0/l;
%rhok=rhok*li;
%phik=phik*li;

phik=ifftshift(rhok);
phi=ifft(phik);
end

