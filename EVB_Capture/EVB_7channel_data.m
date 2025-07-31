clc;clear all;close all;


azi_step=3;
ele_step=3;
AZFOV=360;
ELFOV = 72;
NeleRotation = ELFOV/ele_step+1;
Nazirotation = AZFOV/azi_step+1;

azi_init = -AZFOV/2;
ele_init = -ELFOV;


% Setup Control from Host PC to EVB
evb= serialport('COM10', 115200);
evb.Timeout = 1;


write(evb, 's', 'char');
write(evb,newline, 'char');

GetConfirmFromEVB(evb, 'WaitForChar');
disp('First frame captured')


aziele_cur = [0 0];
aziele_all = [0 0];
RotatorControl_v4(0, ele_init, azi_init);      % 0, el_move, az_move
aziele_cur = aziele_cur+[azi_init ele_init];
aziele_all = [aziele_all;aziele_cur];

iele_start = 1;
iazi_start = 1;


tic
for iele = iele_start:NeleRotation
    if iele==iele_start
        ttt_iazi=iazi_start;
    else
        ttt_iazi=1;
    end

    for iazi=ttt_iazi:Nazirotation
        aziele_cur
        disp(['AZ ' num2str(iazi), '/', num2str(Nazirotation) '    EL ' num2str(iele), '/', num2str(NeleRotation)]);

        % h.AppActivate('putty'); h.SendKeys('a');
        % write(evb, dec2hex(iTx-1), "char");
        % disp(num2str(iTx));
        write(evb, 's', 'char');
        write(evb,newline, 'char');
        
        GetConfirmFromEVB(evb, 'WaitForChar');
        %GetConfirmFromEVB(evb, 'TxReconfigDone');
        %pause(1);
        %Data_ALL(ind_ang,iTx,:,:) = SpectrumAnalyzer(visa,start_freq,stop_freq,n);

        azi = azi_step;
        ele = 0;
        if azi~=0 %Rotate one step in azimuth
            RotatorControl_v4(0, ele, azi);
            aziele_cur = aziele_cur+[azi ele];
            aziele_all = [aziele_all;aziele_cur];
        end
    end

    azi = -iazi*azi_step;
    ele = 0;
    if azi~=0 %Rotate to starting azi to prepare for next ele
        RotatorControl_v4(0, ele, azi);
        aziele_cur = aziele_cur+[azi ele];
        aziele_all = [aziele_all;aziele_cur];
    end

    azi = 0;
    ele = ele_step;
    RotatorControl_v4(0, ele, azi); %rotate one step in elevation
    aziele_cur = aziele_cur+[azi ele];
    aziele_all = [aziele_all;aziele_cur];

end
toc

%Rotate back to original position
azi = -aziele_cur(1);
ele = -aziele_cur(2);
RotatorControl_v4(0, ele, azi);
aziele_cur = aziele_cur+[azi ele];
aziele_all = [aziele_all;aziele_cur];

clear evb;
close all;

