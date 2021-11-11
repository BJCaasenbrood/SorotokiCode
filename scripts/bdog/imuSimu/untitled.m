%clr;
load('thirdSuccess.mat');
D = out.simout.Data;
t = out.simout.time;

p = (D(:,1) - 5)*20;
b = -D(:,3) + 90;

%cla; plot(b)
a1 = [1:371]';
a2 = [372:837]';
a3 = [838:1152]';
a4 = [1153:1622]';
a5 = [1623:1943]';
a6 = [1944:2408]';
a7 = [2409:2730]';
a8 = [2731:3193]';

X = [b(a1);-b(a2)+b(a1(end));
      b(a3);-b(a4)+b(a3(end));
      b(a5);-b(a6)+b(a5(end));
      b(a7);-b(a8)+b(a7(end));];
  
figure;
T = t(1:3193,:);
P = smooth(p(1:3193,:),25);

L = (103-10.5-5)/1e3;
K = smooth((X*pi/180)/L,25); 

plot(P(200:end),K(200:end))
% plot(t(1:3193,:),b1); hold on;
% plot(t(1:3193,:),p(1:3193,:));