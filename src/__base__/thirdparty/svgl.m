function l=svgl(step,p0,p1)
    
    t=(0:step:1)';
    
    l=(1-t).*p0+t.*p1;
end