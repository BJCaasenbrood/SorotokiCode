function y = isomSE3(x)
R = x(1:3,1:3);
p = x(1:3,4);
q = rot2quat(R);

y = [q(:); p(:)];
end