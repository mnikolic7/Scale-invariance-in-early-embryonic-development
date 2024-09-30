function y = curveGauss(p,x)
% function that is used by fitPeaksGauss. 
% curve that is going to be fiteed to the data.
% author: mnikolic@princeton.edu
global Np;
N=Np;
y=zeros(size(x));
for k=1:3:N*3
    x0=p(k);
    wL=p(k+1);
    A=p(k+2);
    y=y+A*gaussian(x,x0,wL);
end


function y = gaussian(x,x0,wL)
    y=1/(wL*sqrt(2*pi))*exp(-0.5*(x-x0).^2./(wL.^2));
end


end