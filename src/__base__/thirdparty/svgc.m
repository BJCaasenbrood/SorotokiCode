function b0123=svgc(step,p0,p1,p2,p3)
    t=(0:step:1)';


    b01=(1-t).*p0+t.*p1;

    b12=(1-t).*p1+t.*p2;

    b012=(1-t).*b01+t.*b12;


    b23=(1-t).*p2+t.*p3;

    b123=(1-t).*b12+t.*b23;


    b0123=(1-t).*b012+t.*b123;
end