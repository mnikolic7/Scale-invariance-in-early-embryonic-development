function [ param,peak_positions,fval ] = fitPeaksGauss( x, y, Npeaks, method, peaks0, varargin )
%This function fits Npeaks Gaussian peaks to the data provided in the
%vector y.
% it takes in the x and y data which have the same dimension (n by 1)
% Npeaks is the number of peaks which you can vary 
% method is the method in which the initial guess will be selected. 
%        default is 'guess' (any string other than the ones below will default to this)
%               for which the program will guess automatically where the peaks are.
%        you can also use 'provide' which will use the parameters for the fit provided 
%               in the next input argument, or
%        'ask' which will ask the user for graphical input for the peaks. 
% peaks0 is the initial guess for the fit parameters, provided that the method is 'provide'
% if you supply any additional input arguments and the last one of them is the string
%   'displayon' in that case this function will also plot the result in the figure 2.
% the output arguments are peak_positions (a Npeaks  by 1 vector of peak positions),
% param (the best fit parameters of the curve, and fval which is the squared error of the fit)
% param is an 3*Npeaks by 1 vector, that contains the center, width, and intensity for each peak 
% in that order: [center1; width1; intenstiy1; center2, width2; intensity2; etc.] 
% This function also handles nans in the x-y data (by ignoring them). (added Jun 01, 2021)
%   author: mnikolic@princeton.edu
%%

fun=@curveGauss;
intensity_constant=1;%/sqrt(2*pi);
global Np;
Np=Npeaks;

plotON=0;
if (nargin>4) && strcmpi(varargin{1},'displayon')
    plotON=1;
end
    
%%
%find peaks for the initial guess of where the centers of the peaks are
[ypeaks, yp_indx]=findpeaks(smooth(y),'NPeaks',Npeaks,'SortStr','descend');


if (length(ypeaks)<Npeaks)
    ypeaks=repmat(max(y),[1,Npeaks]);
    yp_indx=round(linspace(1,length(x),Npeaks));
end

param0=zeros(1,3*Npeaks);
for n=1:Npeaks
    idx_x0=(n-1)*3+1;
    idx_wL=(n-1)*3+2;
    idx_A=(n-1)*3+3;
    
    param0(idx_x0)=x(yp_indx(n));
    param0(idx_wL)=2;
    param0(idx_A)=ypeaks(n)*intensity_constant;
end

%%%%%%%%%%%%%%%%
if strcmpi(method, 'provide')
    param0=peaks0;
%%%%%%%%%%%%%%%%
elseif strcmpi(method,'ask')

    choice='No';
    while (strcmpi(choice, 'No'))
        h=figure('position',[80 80 1000 600]);
        plot(x,y,'o-');
        hold on
        %ask the user to pick the initial parameters graphically (almost
        %never used).
        title(['Please click on top of the ', num2str(Npeaks), ' peaks (left to right)']);
        [x_click,y_click]=ginput(Npeaks);
        
        for kk=1:Npeaks
                idx_x0=(kk-1)*3+1;
                idx_A=(kk-1)*3+3;
                
                param0(idx_x0)=x_click(kk);
                param0(idx_A)=y_click(kk)*intensity_constant;
        end
        
        plot(x, fun(param0,x),'r-','linewidth',2)
        title('Initial fit guess');
        
        choice = questdlg('Is the initial fit good?', ...
            'Good initial fit?','Yes','No','Yes');
        
        close(h);
    end    
end

%clear (potential) nans from the data
x(isnan(y))=[];
y(isnan(y))=[];

dx=x(2)-x(1);
%set initial guess
lb0=[1,        1,   1   ]*dx; %pos, width, intensity
ub0=[max(x),   0.04*max(x),    max(y)*50];
lb=repmat(lb0,1,Npeaks);
ub=repmat(ub0,1,Npeaks);
lb(1:3:end)=param0(1:3:end)-300*dx; %a good liberal estimate of upper and lower bounds that works for pair rule gene curve fitting.
ub(1:3:end)=param0(1:3:end)+300*dx;

options=optimoptions('lsqcurvefit','display','off','algorithm','trust-region-reflective',...
    'MaxFunEvals',1e4, 'maxiter', 1e5, 'TolX', 1e-6, 'TolFun', 1e-6);

[param,fval]=lsqcurvefit(fun,param0,x,y,lb,ub,options);


yfit=curve(param,x);
fval=1./length(x)*sum((y-yfit).^2./(yfit.^2));


if plotON
    figure(2)
    plot(x,y,'bo-');
    hold on

    plot(linspace(x(1),x(end),200),curveGauss(param,linspace(x(1),x(end),200)),'r-');
    for jj=1:Npeaks
        plot(param(1+3*(jj-1))*[ 1 1], [( param(3+3*(jj-1))./(sqrt(2*pi)*param(2+3*(jj-1)) )  ) max(y)],'m*-');
    end
    axis([x(1) x(end) 0 max(y)*1.2]);
    hold off
end
peak_positions=param(1:3:end);
end

