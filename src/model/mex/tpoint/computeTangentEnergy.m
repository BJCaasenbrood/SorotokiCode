function E = computeTangentEnergyFast(gc,PS1)

% compute minimal torus
E = zeros(length(gc),1);

for ii = 1:length(gc)

    g = gc(:,:,ii);
    T = g(1:3,1);
    x = g(1:3,4);
    
    for jj = 1:length(PS1)
        dr = (PS1 - x.');
        R2 = sqrt(diag(dr*dr.')); 
        R2(abs(R2) <= 1e-6) = max(R2);
        
        N = dr./sqrt(sum(dr.^2,2));
        %TangX = cross(repmat(T.',length(R2),1),N);
        TangX(isnan(TangX(:,1)),:) = 0;
        T2 = diag(TangX*TangX.');        
        r = T2./R2;
    end
%     
    E(ii) = sum(r);
end

end