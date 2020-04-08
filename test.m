
t = stepspace(0,0.001,6000);


y = zeros(6000,1);
y2 = zeros(6000,1);

for ii = 1:6000
   y(ii) = st(t(ii),0,2);
   y2(ii) = st(t(ii),3,4);
end

plot(t,y,t,y2)

function y = st(x,a,b)
y = max(0.0,min((x - a) / abs(b - a),1.0));
y =  y * y * (3 - 2 * y);
end