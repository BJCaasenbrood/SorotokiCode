function bus = dev2i2c(device)

switch(device)
    case 'MPRLS',       bus = {'0x18'};
    case 'MPU6050',     bus = {'0x65'};
    case 'LDC1612',     bus = {'0x2A'};
    case 'VEAB',        bus = {'0x62','0x63'};
    otherwise, bus = -1; 
               cout('err',' device not supported by SOROTOKI! \n');
end


end

