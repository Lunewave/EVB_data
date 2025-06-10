clc;clear all;close all;

FOV=10;
azi_step=1;
ele_step=0;

% start_freq = 4.072e9% - halfspan;   
% stop_freq = 4.317e9;% + halfspan; %Calculate frequency span
% n = 21;             %  number of frequency saved


%visalist = visadevlist; %List all connected devices
%ind = find(visalist{:,3} == "Hewlett-Packard"); %Find the device for the Agilent
%visa = visadev(visalist{ind,1}); %Register device in matlab

aziele_cur = [0 0];

% Setup Control from Host PC to EVB
% h= actxserver('WScript.shell');
evb= serialport('COM5', 115200);
evb.Timeout = 1;


RotatorControl_v4(0, 0, -FOV/2);      % 0, el_move, az_move
aziele_cur = [-FOV/2 0];

% write(evb, "s", "char");

tic
for ind_ang=1:(FOV+1)
        ind_ang
           

            % h.AppActivate('putty'); h.SendKeys('a');
            % write(evb, dec2hex(iTx-1), "char");
            % disp(num2str(iTx));
            write(evb, 's', 'char');
            write(evb,newline, 'char');

            GetConfirmFromEVB(evb, 'WaitForChar');
            %GetConfirmFromEVB(evb, 'TxReconfigDone');
            %pause(1);
            %Data_ALL(ind_ang,iTx,:,:) = SpectrumAnalyzer(visa,start_freq,stop_freq,n);


            RotatorControl_v4(0, 0, azi_step);
            aziele_cur(1) = aziele_cur(1) + azi_step;
            aziele_cur(2) = aziele_cur(2) + ele_step;

end
toc

RotatorControl_v4(0, 0, -FOV/2-azi_step);

clear evb;
close all;

