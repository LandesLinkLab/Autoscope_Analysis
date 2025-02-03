function [c,g]=fn_lorentz_fit(x,y,guess,fitnum,lambda_1,lambda_2)
%Lorentz fit
%Fits the spectrum to a sum of one, two, or three Lorentzian peaks
%b1 should be the peak maximum, c1 should be the FWHM

%One Lorentzian
if fitnum ==1
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))',...
         'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a1', 'b1', 'c1'});

%parameters
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 min(x) 0],...
    'Upper',[Inf max(x) Inf],'MaxIter',8000,'tolfun',1e-6,'tolx',1e-6);
st_ = [100 guess 70];  %starting points
set(fo_,'Startpoint',st_);

[c,g] = fit(x,y,f_t,fo_);
end

%%%%% Change the limits if you use these multi-Lorentzian fits in the 
%%%%% future - LSS 021813 

%Two Lorentzian
if fitnum==2
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))+(2*a2/pi).*(c2./(4*(x-b2).^2+c2.^2))');

%parameters
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0 min(x) min(x) 0 0],'Upper',[Inf Inf max(x) max(x) 120 120],'MaxIter',8000);
st_ = [1 1 lambda_1 lambda_2 58.54 62.5];  %starting points
set(fo_,'Startpoint',st_);

[c,g] = fit(x,y,f_t,fo_);
end

%Three Lorentzian
if fitnum ==3
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))+(2*a2/pi).*(c2./(4*(x-b2).^2+c2.^2))+(2*a3/pi).*(c3./(4*(x-b3).^2+c3.^2))');

%parameters
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0 0 400 400 600 0 0 0],'Upper',[Inf Inf Inf 800 800 1200 Inf Inf Inf],'MaxIter',2000);
st_ = [0.5 0.5 0.5 600 650 850 111.6 100 100];  %starting points
set(fo_,'Startpoint',st_);

[c,g] = fit(x,y,f_t,fo_);
end