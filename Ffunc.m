function [F] = Ffunc(h,k,K,,-,-,q_rain)
G = zeros(Nx*Nz,1);
Q = zeros(Nx*Nz,1);
for ii = 1:Nx
    for jj = 1:Nz
        flux_e = ((1-sigma)*k(ii,jj,h(ii,jj))+ sigma*k(ii+1,jj,h(ii+1,jj)))*-(((h(ii+1,jj)+z(jj)) - (h(ii,jj) + z(jj)))*K(ii,jj,east))/delta(ii,jj,east);
        flux_w = -((1-sigma)*k(ii-1,jj,h(ii-1,jj))+ sigma*k(ii,jj,h(ii,jj)))*-(((h(ii,jj)+z(jj)) - (h(ii-1,jj) + z(jj)))*K(ii,jj,west))/delta(ii,jj,west);
        flux_n = ((1-sigma)*k(ii,jj,h(ii,jj))+ sigma*k(ii,jj+1,h(ii,jj+1)))*-(((h(ii,jj+1)+z(jj)) - (h(ii,jj) + z(jj)))*K(ii,jj,north))/delta(ii,jj,north);
        flux_s = -((1-sigma)*k(ii,jj-1,h(ii,jj-1))+ sigma*k(ii,jj,h(ii,jj)))*-(((h(ii,jj)+z(jj)) - (h(ii,jj-1) + z(jj)))*K(ii,jj,south))/delta(ii,jj,south);
        if jj == 1
            flux_s = 0;
        end
        if jj == Nz
            flux_n = q_rain;
        end
        
        if ii==1
            flux_w = 0;
            if z(jj) >= 3 && z(jj) <= 5
                flux_w = -Kc*(Hc - (h(ii,jj) + z(jj)))/xc;
            end
        end
        if ii==Nx
            flux_e = 0;
        end
        G(Nz*(ii-1)+jj) = 1/(Delta(ii,jj,1)*Delta(ii,jj,2))*(flux_e + flux_w + flux_n + flux_s);
        if x(ii)>=0 && x(ii)<=15
            if z(ii) >=8.5 && z(ii)<=10
                Q(Nz*(ii-1)+jj) = -R1(z-8.5)^2/(1.5)^2;
            end
        elseif x(ii) >15 && x(ii)<=100
            if z(ii) >= 7 && z(ii) <=10
                Q(Nz*(ii-1)+jj) = (-0.2*(z(ii)-7)^2)/3^2;
            end
        else
            Q(Nz*(ii-1)+jj) = 0;
        end
    end
end
%hello 


