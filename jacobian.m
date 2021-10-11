function J = jacobian(h_n,F,F_k,options)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% 
% [Nz,Nx] = size(h_n);
% N = Nx*Nz;
% J1 = zeros(N,N);     % Initialise J
% I = eye(N);
% % Variables
% norm_hn = norm(h_n,2);
% 
% if norm_hn == 0
%     epsilon = sqrt(eps);
% else
%     epsilon = sqrt(eps)*norm_hn;
% end
% 
% % Approximate J
% for j = 1:N
%     % determine J at current j
%     J1(:,j) = (F(h_n + epsilon*I(:,j)) - F_k)/epsilon;
% end

%% BANDED
Nx = options.Nx;
Nz = options.Nz;
N = Nx*Nz;
% J = zeros(N,N);     % Initialise J
% I = eye(N);
% Variables
norm_hn = norm(h_n,2);
BW = 2*Nx +1;
s_j = zeros(BW,N);
for ii = 1:BW
    s_j(ii,ii:BW:N) = 1;
end 
%%

J_s_j = zeros(BW,N);
    for j = 1:BW
        % determine epsilon
        if norm_hn == 0
            epsilon = sqrt(eps)/norm(s_j(j,:),2);
        else
            epsilon = (sqrt(eps)*norm_hn)/norm(s_j(j,:),2);
        end
        % determine J_s_j for current j
        J_s_j(j,:) = (F(h_n + epsilon*s_j(j,:)') - F_k)/epsilon;
    end
%%
J = zeros(N,N);
for i = 1:((BW+1)/2)-1
    J([1:i+((BW+1)/2)-1],i) = J_s_j(i,[1:i+((BW+1)/2)-1]);
end

for i = ((BW+1)/2):N - ((BW+1)/2)+1
    if i<=BW
        J([(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW],i) = J_s_j(i,[(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW])';
    else
        J([(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW],i) = J_s_j(mod(i,BW)+(mod(i,BW)==0)*BW,[(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW])';
    end
end

for i = N - ((BW+1)/2)+2:N
    J([(i - (BW+1)/2)+1:N],i) = J_s_j(mod(i,BW)+(mod(i,BW)==0)*BW,[(i - (BW+1)/2)+1:N])';
end

end

