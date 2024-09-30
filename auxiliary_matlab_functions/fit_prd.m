function [peak_position,param,fval] = fit_prd(y, I, pp0)
% fit prd profile with variable number of peaks. for this you must provide
% pp0 [peak positions for each peak you want to fit] 
% depending on the data. pp0 is always provided as relative
% position (units of embryo length). After that this function will
% automatically scale it to max(xdata).
% author: mnikolic@princeton.edu

%% Fit: prd expression profile

ydata=I;
xdata=y';

Npks=length(pp0);
max_X=max(xdata);
max_Y=max(ydata);


%initial position guesses (from mean profile, manually provided).
x00=pp0*max_X;
w00=ones(Npks,1)*0.01;
y00=ones(Npks,1)*0.5*max_Y;
param0=ones(3*Npks,1);


for k=1:Npks
    param0(3*k-2)=x00(k);
    param0(3*k-1)=w00(k);
    param0(3*k)=y00(k);
end


[param,fval]=fitPeaksGauss(xdata,ydata,Npks,'provide',param0,'displayoff');


peak_positions=param(1:3:end);
peak_position=max(peak_positions);

end

