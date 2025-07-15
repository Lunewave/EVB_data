clc;clear all;close all;


for i = 1:100
    % Setup Control from Host PC to EVB
    evb= serialport('COM10', 115200);
    evb.Timeout = 1;
    
    
    write(evb, 's', 'char');
    write(evb,newline, 'char');
    
    GetConfirmFromEVB(evb, 'WaitForChar');
    disp(['Frame ' num2str(i) ' Captured!'])
end

