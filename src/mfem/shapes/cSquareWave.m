
function c = cSquareWave(xc,yc,zc,r,l,N)    
 
if sign(N) == 1
    
c = l/(2*pi*N);
t = 2*pi*N;
R_ = (r/WaveForm(0.5*pi))*0.5;
    
c = @(P) [P(:)*0 + xc, R_*WaveForm(t*P(:)) + yc, c*t*P(:) + zc];

else
    
    c = l/(2*pi*N);
    t = 2*pi*N;
    R_ = (r/WaveForm(0.5*pi))*0.5;
    
    c = @(P) [P(:)*0 + xc, R_*WaveForm(t*P(:)+pi) + yc, -c*t*P(:) + zc];

end

end

function y = WaveForm(x)
delta = 0.02;
y = 2*atan(sin(x(:))/delta)/pi;
end


    