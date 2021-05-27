function h = plotSE3(H)

p = H(1:3,4);
R = H(1:3,1:3);

for ii = 1:3
    plotarrow(p,R(:,ii),ii);
end

end

function plotarrow(p,n,id)
hold on; quiver3(p(1),p(2),p(3), n(1), n(2), n(3),'Color',col(id));
end


