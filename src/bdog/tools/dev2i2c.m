function bus = dev2i2c(device)

switch(device)
    case 'MPRLS',       bus = {'0x18'};
    case 'MPU6050',     bus = {'0x65'};
    case 'LDC1612',     bus = {'0x2A'};
    case 'ADS1013',     bus = {'0x60','0x61'};
    case 'ADS1014',     bus = {'0x60','0x61'};
    case 'ADS1015',     bus = {'0x60','0x61'};
    case 'MP4725',      bus = {'0x48','0x49'};    
        
    otherwise
        bus = -1; 
        cout('err',' device not supported by SOROTOKI! \n');
end

end

