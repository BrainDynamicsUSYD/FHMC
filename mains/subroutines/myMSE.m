function MSE = myMSE(s,dT)
%calculate MSE of a sampler from true mean

%dT = 1e4; %ms
dT=dT+1; %to be consistent with matlab autocorr
%if ~exist('MSE','var')
MSE = 0;
k = 1;
for i=1:1e3:length(s)-dT
    %ss  =  angle(cumsum(  exp(1i*s(i + (1:dT)-1)) ));
    %dss = wrapToPi(ss - true_mean);
    dss = cumsum(s(i + (1:dT)-1))./(1:dT)'; %cum mean
    MSE = MSE + dss.*dss;
    k = k+1;
end
MSE = MSE/k;%(length(s)-dT);
%end
%t = 1:dT;