function y = isomSE3(x)
R = x(1:3,1:3);
p = x(1:3,4);
q = rot2Quat(R);

y = [q(:); p(:)];
end

function q = rot2Quat(R)

    R = squeeze(R);
    T = 1.0 + R(1,1) + R(2,2) + R(3,3);

    if T > 0.00000001

        S  = 0.5 / sqrt(T);
        qw = 0.25 / S;
        qx = (R(3,2) - R(2,3)) * S;
        qy = (R(1,3) - R(3,1)) * S;
        qz = (R(2,1) - R(1,2)) * S;
    elseif R(1,1) > R(2,2) && R(1,1) > R(3,3) 
        S = sqrt( 1.0 + R(1,1) - R(2,2) - R(3,3) ) * 2; 
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S; 
        qz = (R(1,3) + R(3,1)) / S; 
        qw = (R(2,3) - R(3,2)) / S;
    elseif R(2,2) > R(3,3)
        S = sqrt( 1.0 + R(2,2) - R(1,1) - R(3,3) ) * 2; 
        qx = (R(1,2) + R(2,1) ) / S; 
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2) ) / S; 
        qw = (R(1,3) - R(3,1) ) / S;
    else 
        S = sqrt( 1.0 + R(3,3) - R(1,1) - R(2,2) ) * 2; 
        qx = (R(1,3) + R(3,1)) / S; 
        qy = (R(2,3) + R(3,2)) / S; 
        qz = 0.25 * S;
        qw = (R(1,2) - R(2,1)) / S;
    end
    
    q(1) = qw;
    q(2) = qx;
    q(3) = qy;
    q(4) = qz;
  
end
