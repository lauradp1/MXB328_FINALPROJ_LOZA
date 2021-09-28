BW = 43;
N = 20;
s_j = zeros(BW,N);
for ii = 1:BW
    s_j(ii,ii:BW:N) = 1;
end 
%%

J_s_j = zeros(BW,N);
    for j = 1:BW
        % determine epsilon
        if norm_Sn == 0
            epsilon = sqrt(eps)/norm(s_j(j,:),2);
        else
            epsilon = (sqrt(eps)*norm_Sn)/norm(s_j(j,:),2);
        end
        % determine J_s_j for current j
        J_s_j(j,:) = (F(S_n + epsilon*s_j(j,:)') - F(S_n))/epsilon;
    end
%%
J = zeros(N,N);

for i = 1:((BW+1)/2)-1
    J([1:i+((BW+1)/2)-1],i) = J_s_j(i,[1:i+((BW+1)/2)-1]);
end

for i = ((BW+1)/2):N - ((BW+1)/2)+1
    if i<BW
        J([(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW],i) = J_s_j(i,[(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW])';
    else
        J([(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW],i) = J_s_j(mod(i,BW)+1,[(i - (BW+1)/2)+1:(i - (BW+1)/2)+BW])';
    end
end

for i = N - ((BW+1)/2)+2:N
    J([(i - (BW+1)/2)+1:N],i) = J_s_j(mod(i,BW)+1,[(i - (BW+1)/2)+1:N])';
end

