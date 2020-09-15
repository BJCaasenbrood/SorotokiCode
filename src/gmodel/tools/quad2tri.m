function T = quad2tri(E)
T = zeros(2,3);
T(1,:) = [E(1),E(2),E(3)];
T(2,:) = [E(1),E(3),E(4)];
end

