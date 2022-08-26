function [nwse] = NeweyWestDK(e,X,L,constant)
% PURPOSE: computes Newey-West adjusted heteroscedastic-serial
%          consistent standard errors
%
%        Adapted by Guillaume Nolin from the original code by Ian Gow
%---------------------------------------------------
% where: e = T x H matrix of model residuals
%        X = T x k matrix of independant variables
%        L = lag length to use (Default: Newey-West(1994) plug-in
%        procedure)

%        constant = 0: no constant to be added;
%                 = 1: constant term to be added (Default = 1)
%
%        nwse = Newey-West standard errors
%---------------------------------------------------

%% Variables

if nargin < 4 || constant ~= 0
    constant = 1;
end

if ~exist('X','var') && constant ~= 0
    X=ones(size(e,1),1);
end

indexxx = sum(isnan(X),2)==0;
X = X(indexxx,:);
e = e(indexxx,:);

[N,k] = size(X);

if nargin < 3 || L < 0
    % Newey-West (1994) plug-in procedure
    L = floor(4*((N/100)^(2/9)));
end
    
if any(all(X==1,1),2)
    constant=0;
end

if constant == 1
    k = k+1;
    X = [ones(N,1),X];
end

%% Computation
emat=e;
H=size(emat,2);
e=sum(e,2);
evec=reshape(emat',[],1);
Xdiag=kron(eye(H),X);
%Xkeep=X;
%X=repmat(X,[1 H]);

% Q = 0;
% for l = 0:L
%     w_l = 1-l/(L+1);
%     for t = l+1:N
%         if (l==0)   % This calculates the S_0 portion
%             Q = Q  + e(t) ^2 * X(t, :)' * X(t,:);
%         else        % This calculates the off-diagonal terms
%             Q = Q + w_l * e(t) * e(t-l)* ...
%                 (X(t, :)' * X(t-l,:) + X(t-l, :)' * X(t,:));
%         end
%     end
% end
% % Here we average the shocks, then square them, then multiply by the xs
% % (e1+e2)*(x1+x2)=(e1x1+e2x1+e1x2+e2x2)
% % Alternative is to multiply then average
% % (e1x1+e2x2)
% Q = (1/(N-k)) .*Q;

for i=1:size(emat,2)
    h_it_mat(:,:,i)=Xdiag(1+size(X,1)*(i-1):size(X,1)*i,:);
end

for i=1:size(emat,1)
    for j=1:size(emat,2)
        h_it_mat(i,:,j)=emat(i,j)*h_it_mat(i,:,j);
    end
end
h_it_mat_keep=h_it_mat;
h_it_mat=squeeze(sum(h_it_mat,3));

Q2 = 0;
for l = 0:L
    w_l = 1-l/(L+1);
    for t = l+1:N
        if (l==0)   % This calculates the S_0 portion
            Q2 = Q2  + h_it_mat(t,:)'*h_it_mat(t,:);
        else        % This calculates the off-diagonal terms
            Q2 = Q2 + w_l * (h_it_mat(t,:)'*h_it_mat(t-l,:)+h_it_mat(t-l,:)'*h_it_mat(t,:));
        end
    end
end
Q2 = (1/(N-k)) .*Q2;
%Now average Q2;

%nwse = sqrt(diag(N.*((X'*X)\Q/(X'*X))));
nwse = N.*((Xdiag'*Xdiag)\Q2/(Xdiag'*Xdiag));
woo=nwse(k*[1:H],k*[1:H]); % Store the covariance matrix of the stored parameters

end
% Q2 is HK. Average over H then replicate