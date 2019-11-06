function [estimates]= fit_ADC(b_list,vect_SI, start_point)

% --- Create fit "fit 1"
%'Algorithm','Levenberg-Marquardt'
fo = fitoptions('method','NonlinearLeastSquares','Lower',1e-7,'Upper',1);
ok = isfinite(b_list) & isfinite(vect_SI);

set(fo,'Startpoint',start_point);
ft = fittype('exp(-a*x);',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a'});

% Fit this model using new data
cf = fit(b_list(ok)',vect_SI(ok)',ft,fo);
estimates(1)=cf.a;
