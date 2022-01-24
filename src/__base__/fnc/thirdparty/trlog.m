function [w, theta] = trlog(R)

if abs(trace(R) - 3) < 100*eps
    % matrix is identity
    
    w = [0 0 0]';
    theta = 0;
    
elseif abs(trace(R) + 1) < 100*eps
    % tr R = -1
    % rotation by +/- pi, +/- 3pi etc
    
    [mx,k] = max(diag(R));
    I = eye(3,3);
    col = R(:,k) + I(:,k);
    w = col / sqrt(2*(1+mx));
    
    theta = pi;

else
    % general case
    theta = acos( (trace(R)-1)/2 );
    
    skw = (R-R')/2/sin(theta);
    w = vex(skw);   % is a unit vector
    
end
end

