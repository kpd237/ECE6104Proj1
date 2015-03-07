function phi = specpoisson(rho,dx)

eps0=1;
%Take fft and reorder.
rhok=fft(rho);
rhok=fftshift(rhok);

N=length(rhok);
n=linspace(-N/2,N/2-1,N).';

%Calculate Kn^2
kn2=pi*n/N;
Kn2i=(dx/2./sin(kn2*dx)).^2;

phik=rhok./eps0.*Kn2i;%calc charge by dividing by Kn^2
phik=ifftshift(phik)/(N*pi);%put fft back in correct order, scale
phik(1)=0;%2*pi/(dx*N);%eliminate inf from divison by zero.
phi=ifft(phik);
end 
