KX=size(results(1).NL.X,2);
KR=size(results(end).NL.XR,2);
K=KX+KR;
CIparam=norminv(1-(1-CI)/2,0,1);
negshockindex=KX;
posshockindex=KX-1;

for n=1:size(results,2)-1; % loop over LHS variables

for h=0:H
    % heteroscedastic but uncorrelated case
    results(n).NL.variance(:,:,h+1)=var(results(n).NL.res(:,h+1))*inv(results(n).NL.X'*results(n).NL.X);
end

         results(n).NL.variance_beta=squeeze(results(n).NL.variance([KX-1 KX],[KX-1 KX],:));
         results(n).NL.varDK=NeweyWestDK([results(n).NL.res; repmat(results(end).NL.res(:,1),[1 H+1])],blkdiag(results(n).NL.X,results(end).NL.XR),H+1,0);

         %        need to select the right elements from this
         indices=[(posshockindex+K*[0:H]) (negshockindex+K*[0:H]) K*(H+1)-1 K*(H+1)];   

         results(n).NL.varDK_beta=results(n).NL.varDK(indices,indices);

end

for h=0:H
    results(end).NL.variance(:,:,h+1)=var(results(end).NL.res(:,h+1))*inv(results(end).NL.XR'*results(end).NL.XR);
end

results(end).NL.variance_beta=squeeze(results(end).NL.variance([end-1 end],[end-1 end],:));
results(end).NL.varDK=NeweyWestDK(results(end).NL.res,results(end).NL.XR,H+1,0);
results(end).NL.varDK_beta=results(end).NL.varDK(([KR*[1:H+1]-1 KR*[1:H+1]]),([KR*[1:H+1]-1 KR*[1:H+1]]));

% delta method 
for n=1:size(results,2)-1; % loop over LHS variables

    % scaled impulse responses
    results(n).NL.scaled.delta_pos=[eye(H+1)/results(end).NL.betafull(1,end-1) zeros(H+1) -results(n).NL.betafull(:,end-1)/(results(end).NL.betafull(1,end-1)^2) zeros(H+1,1)];
    results(n).NL.scaled.delta_neg=[zeros(H+1) eye(H+1) zeros(H+1,2)];  
    results(n).NL.scaled.delta_diff=results(n).NL.scaled.delta_pos-results(n).NL.scaled.delta_neg;

    results(n).NL.scaled.var_pos=results(n).NL.scaled.delta_pos*results(n).NL.varDK_beta*results(n).NL.scaled.delta_pos';
    results(n).NL.scaled.var_neg=results(n).NL.scaled.delta_neg*results(n).NL.varDK_beta*results(n).NL.scaled.delta_neg';
    results(n).NL.scaled.var_diff=results(n).NL.scaled.var_pos-results(n).NL.scaled.var_neg;;

    results(n).NL.scaled.SE_pos=sqrt(diag(results(n).NL.scaled.var_pos));
    results(n).NL.scaled.SE_neg=sqrt(diag(results(n).NL.scaled.var_neg));
    results(n).NL.scaled.SE_diff=sqrt(diag(results(n).NL.scaled.var_diff));

    results(n).NL.scaled.tstats_pos=results(n).NL.scaled.coef(:,1)./results(n).NL.scaled.SE_pos; % T stats defined as standard
    results(n).NL.scaled.tstats_neg=results(n).NL.scaled.coef(:,2)./results(n).NL.scaled.SE_neg; % T stats defined as standard
    results(n).NL.scaled.tstats_diff=results(n).NL.scaled.coef_diff./sqrt(diag(results(n).NL.scaled.var_diff));

    results(n).NL.scaled.CIhigh_pos =  results(n).NL.scaled.coef(:,1) + CIparam*results(n).NL.scaled.SE_pos; % 90% confidence intervals
    results(n).NL.scaled.CIhigh_neg =  results(n).NL.scaled.coef(:,2) + CIparam*results(n).NL.scaled.SE_neg; % 90% confidence intervals
    results(n).NL.scaled.CIlow_pos =  results(n).NL.scaled.coef(:,1) - CIparam*results(n).NL.scaled.SE_pos; % 90% confidence intervals
    results(n).NL.scaled.CIlow_neg =  results(n).NL.scaled.coef(:,2) - CIparam*results(n).NL.scaled.SE_neg; % 90% confidence intervals

    % smoothed impulse responses
    results(n).NL.smoothed.delta_pos=[weight/results(end).NL.betafull(1,end-1) zeros(H+1) -weight*results(n).NL.betafull(:,end-1)/(results(end).NL.betafull(1,end-1)^2) zeros(H+1,1)];
    results(n).NL.smoothed.delta_neg=[zeros(H+1) weight zeros(H+1,2)];  
    results(n).NL.smoothed.delta_diff=results(n).NL.smoothed.delta_pos-results(n).NL.smoothed.delta_neg;

    results(n).NL.smoothed.var_pos=results(n).NL.smoothed.delta_pos*results(n).NL.varDK_beta*results(n).NL.smoothed.delta_pos';
    results(n).NL.smoothed.var_neg=results(n).NL.smoothed.delta_neg*results(n).NL.varDK_beta*results(n).NL.smoothed.delta_neg';
    results(n).NL.smoothed.var_diff=results(n).NL.smoothed.delta_diff*results(n).NL.varDK_beta*results(n).NL.smoothed.delta_diff';

    results(n).NL.smoothed.SE_pos=sqrt(diag(results(n).NL.smoothed.var_pos));
    results(n).NL.smoothed.SE_neg=sqrt(diag(results(n).NL.smoothed.var_neg));
    results(n).NL.smoothed.SE_diff=sqrt(diag(results(n).NL.smoothed.var_diff));

    results(n).NL.smoothed.tstats_pos=results(n).NL.smoothed.coef(:,1)./results(n).NL.smoothed.SE_pos; % T stats defined as standard
    results(n).NL.smoothed.tstats_neg=results(n).NL.smoothed.coef(:,2)./results(n).NL.smoothed.SE_neg; % T stats defined as standard
    results(n).NL.smoothed.tstats_diff=results(n).NL.smoothed.coef_diff./sqrt(diag(results(n).NL.smoothed.var_diff));

    results(n).NL.smoothed.CIhigh_pos =  results(n).NL.smoothed.coef(:,1) + CIparam*results(n).NL.smoothed.SE_pos; % 90% confidence intervals
    results(n).NL.smoothed.CIhigh_neg =  results(n).NL.smoothed.coef(:,2) + CIparam*results(n).NL.smoothed.SE_neg; % 90% confidence intervals
    results(n).NL.smoothed.CIlow_pos =  results(n).NL.smoothed.coef(:,1) - CIparam*results(n).NL.smoothed.SE_pos; % 90% confidence intervals
    results(n).NL.smoothed.CIlow_neg =  results(n).NL.smoothed.coef(:,2) - CIparam*results(n).NL.smoothed.SE_neg; % 90% confidence intervals

    % cumulated impulse responses
    results(n).NL.cum.delta_pos=[cumweight/results(end).NL.betafull(1,end-1) zeros(H+1) -cumweight*results(n).NL.betafull(:,end-1)/(results(end).NL.betafull(1,end-1)^2) zeros(H+1,1)];
    results(n).NL.cum.delta_neg=[zeros(H+1) cumweight zeros(H+1,2)];  
    results(n).NL.cum.delta_diff=results(n).NL.cum.delta_pos-results(n).NL.cum.delta_neg;

    results(n).NL.cum.var_pos=results(n).NL.cum.delta_pos*results(n).NL.varDK_beta*results(n).NL.cum.delta_pos';
    results(n).NL.cum.var_neg=results(n).NL.cum.delta_neg*results(n).NL.varDK_beta*results(n).NL.cum.delta_neg';
    results(n).NL.cum.var_diff=results(n).NL.cum.delta_diff*results(n).NL.varDK_beta*results(n).NL.cum.delta_diff';

    results(n).NL.cum.SE_pos=sqrt(diag(results(n).NL.cum.var_pos));
    results(n).NL.cum.SE_neg=sqrt(diag(results(n).NL.cum.var_neg));
    results(n).NL.cum.SE_diff=sqrt(diag(results(n).NL.cum.var_diff));

    results(n).NL.cum.tstats_pos=results(n).NL.cum.coef(:,1)./results(n).NL.cum.SE_pos; % T stats defined as standard
    results(n).NL.cum.tstats_neg=results(n).NL.cum.coef(:,2)./results(n).NL.cum.SE_neg; % T stats defined as standard
    results(n).NL.cum.tstats_diff=results(n).NL.cum.coef_diff./sqrt(diag(results(n).NL.cum.var_diff));

    results(n).NL.cum.CIhigh_pos =  results(n).NL.cum.coef(:,1) + CIparam*results(n).NL.cum.SE_pos; % 90% confidence intervals
    results(n).NL.cum.CIhigh_neg =  results(n).NL.cum.coef(:,2) + CIparam*results(n).NL.cum.SE_neg; % 90% confidence intervals
    results(n).NL.cum.CIlow_pos =  results(n).NL.cum.coef(:,1) - CIparam*results(n).NL.cum.SE_pos; % 90% confidence intervals
    results(n).NL.cum.CIlow_neg =  results(n).NL.cum.coef(:,2) - CIparam*results(n).NL.cum.SE_neg; % 90% confidence intervals
% end

end

% policy variable
% scaled impulse responses
results(end).NL.scaled.delta_pos=[eye(H+1)/results(end).NL.betafull(1,end-1) zeros(H+1)];
results(end).NL.scaled.delta_pos(:,1)=results(end).NL.scaled.delta_pos(:,1)-results(end).NL.betafull(:,end-1)/(results(end).NL.betafull(1,end-1)^2);
results(end).NL.scaled.delta_neg=[zeros(H+1) eye(H+1)];
results(end).NL.scaled.delta_diff=results(end).NL.scaled.delta_pos-results(end).NL.scaled.delta_neg;

results(end).NL.scaled.var_pos=results(end).NL.scaled.delta_pos*results(end).NL.varDK_beta*results(end).NL.scaled.delta_pos';
results(end).NL.scaled.var_neg=results(end).NL.scaled.delta_neg*results(end).NL.varDK_beta*results(end).NL.scaled.delta_neg';
results(end).NL.scaled.var_diff=results(end).NL.scaled.delta_diff*results(end).NL.varDK_beta*results(end).NL.scaled.delta_diff';

results(end).NL.scaled.SE_pos=sqrt(diag(results(end).NL.scaled.var_pos));
results(end).NL.scaled.SE_neg=sqrt(diag(results(end).NL.scaled.var_neg));
results(end).NL.scaled.SE_diff=sqrt(diag(results(end).NL.scaled.var_diff));

results(end).NL.scaled.tstats_pos=results(end).NL.scaled.coef(:,1)./results(end).NL.scaled.SE_pos; % T stats defined as standard
results(end).NL.scaled.tstats_neg=results(end).NL.scaled.coef(:,2)./results(end).NL.scaled.SE_neg; % T stats defined as standard
results(end).NL.scaled.tstats_diff=results(end).NL.scaled.coef_diff./sqrt(diag(results(end).NL.scaled.var_diff));

results(end).NL.scaled.CIhigh_pos =  results(end).NL.scaled.coef(:,1) + CIparam*results(end).NL.scaled.SE_pos; % 90% confidence intervals
results(end).NL.scaled.CIhigh_neg =  results(end).NL.scaled.coef(:,2) + CIparam*results(end).NL.scaled.SE_neg; % 90% confidence intervals
results(end).NL.scaled.CIlow_pos =  results(end).NL.scaled.coef(:,1) - CIparam*results(end).NL.scaled.SE_pos; % 90% confidence intervals
results(end).NL.scaled.CIlow_neg =  results(end).NL.scaled.coef(:,2) - CIparam*results(end).NL.scaled.SE_neg; % 90% confidence intervals

% smoothed impulse responses
results(end).NL.smoothed.delta_pos=[weight/results(end).NL.betafull(1,end-1) zeros(H+1)];
results(end).NL.smoothed.delta_pos(:,1)=results(end).NL.smoothed.delta_pos(:,1)-weight*results(end).NL.betafull(:,end-1)/(results(end).NL.betafull(1,end-1)^2);
results(end).NL.smoothed.delta_neg=[zeros(H+1) weight];
results(end).NL.smoothed.delta_diff=results(end).NL.smoothed.delta_pos-results(end).NL.smoothed.delta_neg;

results(end).NL.smoothed.var_pos=results(end).NL.smoothed.delta_pos*results(end).NL.varDK_beta*results(end).NL.smoothed.delta_pos';
results(end).NL.smoothed.var_neg=results(end).NL.smoothed.delta_neg*results(end).NL.varDK_beta*results(end).NL.smoothed.delta_neg';
results(end).NL.smoothed.var_diff=results(end).NL.smoothed.delta_diff*results(end).NL.varDK_beta*results(end).NL.smoothed.delta_diff';

results(end).NL.smoothed.SE_pos=sqrt(diag(results(end).NL.smoothed.var_pos));
results(end).NL.smoothed.SE_neg=sqrt(diag(results(end).NL.smoothed.var_neg));
results(end).NL.smoothed.SE_diff=sqrt(diag(results(end).NL.smoothed.var_diff));

results(end).NL.smoothed.tstats_pos=results(end).NL.smoothed.coef(:,1)./results(end).NL.smoothed.SE_pos; % T stats defined as standard
results(end).NL.smoothed.tstats_neg=results(end).NL.smoothed.coef(:,2)./results(end).NL.smoothed.SE_neg; % T stats defined as standard
results(end).NL.smoothed.tstats_diff=results(end).NL.smoothed.coef_diff./sqrt(diag(results(end).NL.smoothed.var_diff));

results(end).NL.smoothed.CIhigh_pos =  results(end).NL.smoothed.coef(:,1) + CIparam*results(end).NL.smoothed.SE_pos; % 90% confidence intervals
results(end).NL.smoothed.CIhigh_neg =  results(end).NL.smoothed.coef(:,2) + CIparam*results(end).NL.smoothed.SE_neg; % 90% confidence intervals
results(end).NL.smoothed.CIlow_pos =  results(end).NL.smoothed.coef(:,1) - CIparam*results(end).NL.smoothed.SE_pos; % 90% confidence intervals
results(end).NL.smoothed.CIlow_neg =  results(end).NL.smoothed.coef(:,2) - CIparam*results(end).NL.smoothed.SE_neg; % 90% confidence intervals

% cumulated impulse responses
results(end).NL.cum.delta_pos=[cumweight/results(end).NL.betafull(1,end-1) zeros(H+1)];
results(end).NL.cum.delta_pos(:,1)=results(end).NL.smoothed.delta_pos(:,1)-cumweight*results(end).NL.betafull(:,end-1)/(results(end).NL.betafull(1,end-1)^2);
results(end).NL.cum.delta_neg=[zeros(H+1) cumweight];
results(end).NL.cum.delta_diff=results(end).NL.cum.delta_pos-results(end).NL.cum.delta_neg;

results(end).NL.cum.var_pos=results(end).NL.cum.delta_pos*results(end).NL.varDK_beta*results(end).NL.cum.delta_pos';
results(end).NL.cum.var_neg=results(end).NL.cum.delta_neg*results(end).NL.varDK_beta*results(end).NL.cum.delta_neg';
results(end).NL.cum.var_diff=results(end).NL.cum.delta_diff*results(end).NL.varDK_beta*results(end).NL.cum.delta_diff';

results(end).NL.cum.SE_pos=sqrt(diag(results(end).NL.cum.var_pos));
results(end).NL.cum.SE_neg=sqrt(diag(results(end).NL.cum.var_neg));
results(end).NL.cum.SE_diff=sqrt(diag(results(end).NL.cum.var_diff));

results(end).NL.cum.tstats_pos=results(end).NL.cum.coef(:,1)./results(end).NL.cum.SE_pos; % T stats defined as standard
results(end).NL.cum.tstats_neg=results(end).NL.cum.coef(:,2)./results(end).NL.cum.SE_neg; % T stats defined as standard
results(end).NL.cum.tstats_diff=results(end).NL.cum.coef_diff./sqrt(diag(results(end).NL.cum.var_diff));

results(end).NL.cum.CIhigh_pos =  results(end).NL.cum.coef(:,1) + CIparam*results(end).NL.cum.SE_pos; % 90% confidence intervals
results(end).NL.cum.CIhigh_neg =  results(end).NL.cum.coef(:,2) + CIparam*results(end).NL.cum.SE_neg; % 90% confidence intervals
results(end).NL.cum.CIlow_pos =  results(end).NL.cum.coef(:,1) - CIparam*results(end).NL.cum.SE_pos; % 90% confidence intervals
results(end).NL.cum.CIlow_neg =  results(end).NL.cum.coef(:,2) - CIparam*results(end).NL.cum.SE_neg; % 90% confidence intervals
