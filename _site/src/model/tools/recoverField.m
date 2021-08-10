function [u,X] = recoverField(Model,k)
P = Model.get('Phi');
X = linspace(0,Model.get('Sdomain'),250);
Q = Model.q;
u = [];

for ii = 1:length(X)
    Phi = P(X(ii));
    u = vappend(u,(Phi*Q(k,:).').');
end

end

