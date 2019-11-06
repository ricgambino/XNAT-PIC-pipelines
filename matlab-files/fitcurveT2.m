 function [estimates, model,exitflag] = fitcurveT2(vect_echoes, vect_analysis, start_point)

% giugno 2009
%modificato con metodo più robusto per il calcolo della bontà del fitting

% Call fminsearch with a random starting point.

model = @T2_expfun;
[estimates,fval,exitflag] = fminsearch(model, start_point);

%[estimates,fval,exitflag] = fminsearch(model, start_point, optimset('Display','final', 'TolX', 1e-8, 'MaxIter', 1000, 'TolFun', 1e-7) );

% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
 

function [sse, FittedCurve] = T2_expfun(params)
        M0 = params(1);
        T2 = params(2);
        off = params(3);
        %Mxy(t)=Mxy(0)(exp(-t/T2))
        FittedCurve = M0.* exp(-(vect_echoes)/(T2))+off;
        %ErrorVector = FittedCurve - vect_analysis;
        %sse = sum(ErrorVector.^2);
            ss_reg=(vect_analysis(:)-FittedCurve(:)).^2;
            ss_tot=(vect_analysis(:)-mean(vect_analysis(:))).^2;
            SS_reg=sum(ss_reg(:));
            SS_tot=sum(ss_tot(:));
            sse=(SS_reg/SS_tot);
        
        
end
 end