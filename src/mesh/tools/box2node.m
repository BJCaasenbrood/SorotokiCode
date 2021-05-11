function Nodes = box2node(B)

if numel(B) == 4
    Nodes = [B(1),B(3);
        B(2),B(3);
        B(2),B(4);
        B(1),B(4)];
else
    Nodes = [B(1),B(3),B(5);
        B(1),B(3),B(6);
        B(1),B(4),B(5);
        B(1),B(4),B(6);
        B(2),B(3),B(5);
        B(2),B(3),B(6);
        B(2),B(4),B(5);
        B(2),B(4),B(6)];
end
end