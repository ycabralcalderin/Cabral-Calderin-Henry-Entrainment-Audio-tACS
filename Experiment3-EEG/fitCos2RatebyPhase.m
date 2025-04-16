function    [fit, RESNORM,RESIDUAL,prefPhase]=fitCos2RatebyPhase (C,X,params)

%parameters for fitting (start always with middle values)

%Fitting
[fit, RESNORM,RESIDUAL] = lsqcurvefit(@cosfit_rad,params,[X X]',C,[0 0 0],[2*pi 1 1]);
Xnew = 0:2*pi/1000:2*pi-2*pi/1000; %this is my phase vector
y = fit(2) + fit(3).*(cos(Xnew + fit(1))); % if you want to plot the predicted function
[~,pos]=max(y);
prefPhase=Xnew(pos);



end