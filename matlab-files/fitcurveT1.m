 function [estimates, model,exitflag] = fitcurveT1(vect_echoes, vect_analysis, start_point)

 
% Call fminsearch with a random starting point.

model = @T1_expfun;
[estimates,fval,exitflag] = fminsearch(model, start_point);

% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
 

function [sse, FittedCurve] = T1_expfun(params)
        MZ = params(1);
        T1 = params(2);
        off = params(3);
        %Mz(t) = Mzeq(1-exp(-t/T1))
        FittedCurve = (MZ.* (1- exp(-(vect_echoes)/(T1))))+off;
        ErrorVector = FittedCurve - vect_analysis;
        sse = sum(ErrorVector .^ 2);
        
end
 end