function RGB = Normal2RGB(N,Diffuse)

RGB(:,1) = 0.5 + 0.5*(N(:,1));
RGB(:,2) = 0.5 + 0.5*(N(:,2));
RGB(:,3) = 0.5 + 0.5*softabs(N(:,3));

end