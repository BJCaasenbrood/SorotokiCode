function ver = BuildVersion

Inception = '30-feb-18';
NumDays = datenum(date) - datenum(Inception);

c1 = 200; c2 = 30; c3 = 5;

a = mod(NumDays,c1);
a_ = ceil((NumDays - a)/c1);
b = mod(a,c2);
b_ = ceil((a - b)/c2);
c = mod(b,c3);
c_ = ceil((b - c)/c3);

ver = [num2str(a_),'.0',num2str(b_),'.0',num2str(c_)];
end

