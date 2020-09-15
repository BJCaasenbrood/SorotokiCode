function [A,J1,J2,J3]= JRectangle(b,h)
A = b*h;
J1 = (1/12)*b*h*(b^2 + h^2);
J2 = (1/12)*b*h^3;
J3 = (1/12)*h*b^3;
end

%https://www.efunda.com/math/areas/rectangle.cfm