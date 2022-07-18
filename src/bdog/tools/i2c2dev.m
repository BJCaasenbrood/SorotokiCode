function dev = i2c2dev(address)

switch(address)
    case 18,            dev = 'MPRLS';
    case '0x18',        dev = 'MPRLS';
    case 65,            dev = 'MPU6050';
    case '0x65',        dev = 'MPU6050';  
    case '0x2A',        dev = 'LDC1612';  
    case 60,            dev = 'ADS10xx';     
    case 61,            dev = 'ADS10xx';         
    case '0x60',        dev = 'ADS10xx';  
    case '0x61',        dev = 'ADS10xx';  
    case 48,            dev = 'MP4725';     
    case 49,            dev = 'MP4725';         
    case '0x48',        dev = 'MP4725';  
    case '0x49',        dev = 'MP4725';          

    otherwise
       dev = []; 
       cout('err',' device not supported by SOROTOKI! \n');
       
end

end

