%% block bootstrap
blocks=ceil(samplelength/blocklength);
    Rsample=R(startS:endS);

for n=1:size(LHS,2)

    results(n).lin.betaBS=[];
    results(n).lin.betaBSR=[];
    if runstate==1
        results(n).NL.betaBS=zeros(2,21,B);
        results(n).NL.betaBSR=zeros(2,1,B);
    end

    LHSmat=[];
    for h=0:H
        LHSmat=[LHSmat LHS(startS+h:endS+h,n)];
    end

    % draw correlated samples
    for b=1:B
        % pick int(T/L) indices between 1 and L observations from the end of the sample 
        % concatenate the subsamples L periods long starting at each
        % index
        blockindices=randi(samplelength-blocklength+1,blocks,1);
        LHSmatb=[];
        Xb=[];
        XRb=[];
        Rb=[];

        for j=1:blocks
           LHSmatb=[LHSmatb; LHSmat(blockindices(j):blockindices(j)+blocklength-1,:)];
           Xb=[Xb; results(n).lin.X(blockindices(j):blockindices(j)+blocklength-1,:)];
           XRb=[XRb; results(size(LHS,2)+1).lin.XR(blockindices(j):blockindices(j)+blocklength-1,:)];
           Rb=[Rb; Rsample(blockindices(j):blockindices(j)+blocklength-1,:)];
        end

        betaBS=inv(Xb'*Xb)*Xb'*LHSmatb;
        betaBSR=inv(XRb'*XRb)*XRb'*Rb;
        results(n).lin.betaBS=[results(n).lin.betaBS; betaBS(end,:)];
        results(n).lin.betaBSR=[results(n).lin.betaBSR; betaBSR(end,1)];
        
        if runstate==1
            XNLb=[];
            XRNLb=[];
            for j=1:blocks
               XNLb=[XNLb; results(n).NL.XNL(blockindices(j):blockindices(j)+blocklength-1,:)];
               XRNLb=[XRNLb; results(size(LHS,2)+1).NL.XRNL(blockindices(j):blockindices(j)+blocklength-1,:)];
            end
            betaBSNL=inv(XNLb'*XNLb)*XNLb'*LHSmatb;
            betaBSRNL=inv(XRNLb'*XRNLb)*XRNLb'*Rb;

            results(n).NL.betaBS(:,:,b)=betaBSNL([end-(end-includetrend)/2 end],:);
            results(n).NL.betaBSR(:,:,b)=betaBSRNL([end-(end-includetrend)/2 end],:);

        end

    end         % close loop over bootstrap iterations


    % transform to scaled, smoothed and cumulative parameters
    results(n).lin.scaled.coefBS=results(n).lin.betaBS./repmat(results(n).lin.betaBSR,[1 21]); % Scale by first period impact on Fed Funds
    results(n).lin.scaled.coefBS_mean=squeeze(mean(results(n).lin.scaled.coefBS));

    results(n).lin.smoothed.coefBS=results(n).lin.scaled.coefBS*weight';
    results(n).lin.smoothed.coefBS_mean=squeeze(mean(results(n).lin.smoothed.coefBS));

    results(n).lin.cum.coefBS=results(n).lin.scaled.coefBS*weight';
    results(n).lin.cum.coefBS_mean=squeeze(mean(results(n).lin.cum.coefBS));


    if runstate==1
        results(n).NL.scaled.coefBS=results(n).NL.betaBS./repmat(results(n).NL.betaBSR,[1 21 1]); % Scale by first period impact on Fed Funds
        results(n).NL.scaled.coefBS_mean=mean(results(n).NL.scaled.coefBS,3);
                
        results(n).NL.smoothed.coefBS(1,:,:)=weight*squeeze(results(n).NL.scaled.coefBS(1,:,:));
        results(n).NL.smoothed.coefBS(2,:,:)=weight*squeeze(results(n).NL.scaled.coefBS(2,:,:));
        results(n).NL.smoothed.coefBS_mean=mean(results(n).NL.smoothed.coefBS,3);
        
        
        results(n).NL.cum.coefBS(1,:,:)=cumweight*squeeze(results(n).NL.scaled.coefBS(1,:,:));
        results(n).NL.cum.coefBS(2,:,:)=cumweight*squeeze(results(n).NL.scaled.coefBS(2,:,:));
        results(n).NL.cum.coefBS_mean=mean(results(n).NL.cum.coefBS,3);

        results(n).NL.scaled.betadiff=(results(n).NL.scaled.coefBS(1,:,:)<=results(n).NL.scaled.coefBS(2,:,:));
        results(n).NL.scaled.betadiffshare=mean(results(n).NL.scaled.betadiff,3);
        results(n).NL.scaled.betadifft=-norminv(results(n).NL.scaled.betadiffshare,0,1);

        results(n).NL.smoothed.betadiff=(results(n).NL.smoothed.coefBS(1,:,:)<=results(n).NL.smoothed.coefBS(2,:,:));
        results(n).NL.smoothed.betadiffshare=mean(results(n).NL.smoothed.betadiff,3);
        results(n).NL.smoothed.betadifft=-norminv(results(n).NL.smoothed.betadiffshare,0,1);

        results(n).NL.cum.betadiff=(results(n).NL.cum.coefBS(1,:,:)<=results(n).NL.cum.coefBS(2,:,:));
        results(n).NL.cum.betadiffshare=mean(results(n).NL.cum.betadiff,3);
        results(n).NL.cum.betadifft=-norminv(results(n).NL.cum.betadiffshare,0,1);

    end
end

% now do policy variable
n=n+1;
results(n).lin.betaBS=[];
if runstate==1
    results(n).NL.betaBS=zeros(2,21,B);
end
Rsample=R(startS:endS);
Rmat=[];
for h=0:H
    Rmat=[Rmat R(startS+h:endS+h)];
end

results(n).lin.betaBS=[];
% draw correlated samples
for b=1:B
    blockindices=randi(samplelength-blocklength+1,blocks,1);
    Rmatb=[];
    Rb=[];
    XRb=[];

    for j=1:blocks
       Rmatb=[Rmatb; Rmat(blockindices(j):blockindices(j)+blocklength-1,:)];
       XRb=[XRb; results(size(LHS,2)+1).lin.XR(blockindices(j):blockindices(j)+blocklength-1,:)];
    end

    betaBSR=inv(XRb'*XRb)*XRb'*Rmatb;
    results(n).lin.betaBS=[results(n).lin.betaBS; betaBSR(end,:)];

    if runstate==1
        XRNLb=[];
        for j=1:blocks
           XRNLb=[XRNLb; results(size(LHS,2)+1).NL.XRNL(blockindices(j):blockindices(j)+blocklength-1,:)];
        end
        betaBSRNL=inv(XRNLb'*XRNLb)*XRNLb'*Rmatb;
        results(n).NL.betaBS(:,:,b)=betaBSRNL([end-(end-includetrend)/2 end],:);

    end

end         % close loop over bootstrap iterations


% transform to scaled, smoothed and cumulative parameters
results(n).lin.scaled.coefBS=results(n).lin.betaBS./repmat(results(n).lin.betaBS(:,1),[1 21]); % Scale by first period impact on Fed Funds
results(n).lin.scaled.coefBS_mean=squeeze(mean(results(n).lin.scaled.coefBS));

results(n).lin.smoothed.coefBS=results(n).lin.scaled.coefBS*weight';
results(n).lin.smoothed.coefBS_mean=squeeze(mean(results(n).lin.smoothed.coefBS));

results(n).lin.cum.coefBS=results(n).lin.scaled.coefBS*cumweight';
results(n).lin.cum.coefBS_mean=squeeze(mean(results(n).lin.cum.coefBS));

if runstate==1
    results(n).NL.scaled.coefBS=results(n).NL.betaBS./repmat(results(n).NL.betaBS(:,1,:),[1 21 1]); % Scale by first period impact on Fed Funds
    results(n).NL.scaled.coefBS_mean=mean(results(n).NL.scaled.coefBS,3);

    results(n).NL.smoothed.coefBS(1,:,:)=weight*squeeze(results(n).NL.scaled.coefBS(1,:,:));
    results(n).NL.smoothed.coefBS(2,:,:)=weight*squeeze(results(n).NL.scaled.coefBS(2,:,:));
    results(n).NL.smoothed.coefBS_mean=mean(results(n).NL.smoothed.coefBS,3);

    results(n).NL.cum.coefBS(1,:,:)=cumweight*squeeze(results(n).NL.scaled.coefBS(1,:,:));
    results(n).NL.cum.coefBS(2,:,:)=cumweight*squeeze(results(n).NL.scaled.coefBS(2,:,:));
    results(n).NL.cum.coefBS_mean=mean(results(n).NL.cum.coefBS,3);

    results(n).NL.scaled.betadiff=(results(n).NL.scaled.coefBS(1,:,:)<=results(n).NL.scaled.coefBS(2,:,:));
    results(n).NL.scaled.betadiffshare=mean(results(n).NL.scaled.betadiff,3);
    results(n).NL.scaled.betadifft=-norminv(results(n).NL.scaled.betadiffshare,0,1);

    results(n).NL.smoothed.betadiff=(results(n).NL.smoothed.coefBS(1,:,:)<=results(n).NL.smoothed.coefBS(2,:,:));
    results(n).NL.smoothed.betadiffshare=mean(results(n).NL.smoothed.betadiff,3);
    results(n).NL.smoothed.betadifft=-norminv(results(n).NL.smoothed.betadiffshare,0,1);

    results(n).NL.cum.betadiff=(results(n).NL.cum.coefBS(1,:,:)<=results(n).NL.cum.coefBS(2,:,:));
    results(n).NL.cum.betadiffshare=mean(results(n).NL.cum.betadiff,3);
    results(n).NL.cum.betadifft=-norminv(results(n).NL.cum.betadiffshare,0,1);

end

clear results.NL.smoothed.coefBS;
clear results.NL.scaled.coefBS;
clear results.NL.cum.coefBS;
