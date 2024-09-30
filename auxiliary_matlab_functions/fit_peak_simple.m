function [peak_position,param] = fit_peak_simple(x, I)
% fit single peak.
% author: mnikolic@princeton.edu


%% Fit: 

ydata=I;
xdata=x';

max_X=max(xdata);
min_X=min(xdata);
max_Y=max(ydata);
min_Y=min(ydata);

[~,loc]=max(smooth(ydata,200));
center=xdata(loc(1));
param0=[center, 0.005*max_X 1*max_Y, 0];
fun=@gaussian_with_bak;
lb=[xdata(1),        0.0001*(max_X-min_X),   0.01*(max_Y-min_Y), -max_Y]; %pos, width, intensity, bak
ub=[max_X,   1*(max_X-min_X),    max_Y*50, max_Y];

options=optimoptions('lsqcurvefit','display','off','algorithm','trust-region-reflective',...
    'MaxFunEvals',1e4, 'maxiter', 1e5, 'TolX', 1e-6, 'TolFun', 1e-6);


[param_final,fval]=lsqcurvefit(fun,param0,xdata,ydata,lb,ub,options);



peak_position=param_final(1);
param=param_final;
end

function y = gaussian_with_bak(p,x)
    wL=p(2);
    I=p(3);
    x0=p(1);
    B=p(4);
    y=I/(wL*sqrt(2*pi))*exp(-0.5*(x-x0).^2./(wL.^2))+B;
end

