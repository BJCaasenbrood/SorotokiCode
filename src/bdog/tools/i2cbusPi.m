function [SCL,SDA] = i2cbusPi(bus)

switch(bus)
    case 1, SCL = 3; SDA = 2;
    case 0, SCL = 1; SDA = 0;
    case 3, SCL = 5; SDA = 4;
    case 4, SCL= 9; SDA = 8;
    case 5, SCL = 13; SDA = 12;
    case 6, SCL = 23; SDA = 22;
    otherwise, cout('err',' i2c bus not defined! \n');
end

end
