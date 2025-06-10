function RotatorControl_v4(Hdegree,Vdegree,Rdegree)
%% Rotator Chontrol
% Input: degree values to rotate; 
% Usage Example: RotatorControl(30,60,90);Rotate 30 degrees Horizonally and
% 60 vertically and Roll 90 degrees
% Note: It is necessary to change the serial port number COM. 
%=========================%
% Ning Cao-Oct2019-Tucson
% NC 20200902:
%=========================%
%% Note; For V3 Rotator
% the ONLY thing you really need to change is the gear ratio to 72:1
% instead of the 14.4:1 for V2
% Also , for the current settings for both motors, set the value to 200,75,10    
% just be sure not to use the old turntable with those current settings , 
% it will make the motors really hot and possibly burn them up.
% Hdegree= 10;
% Vdegree=10;
% Rdegree = 10;
%%
% The resolution of Z stage is 90:1 , the other stages (XY) are 72:1 ,  
% I think if you change your Z formula to 16*10*5 you should get 800 …..   
% …. That stage is a NATIVE  50 steps per degree so with 1/16th microstepping
% you get 800 steps per degree.
%% Try TM toolbox
% In Command Window: tmtool
%% COM number
COMNum = 'COM4';
%% In order to avoid damage the Rotator, the vertical rotation degrees are limited to be
VdegreeLimit = [-30 30];
if Vdegree<VdegreeLimit(1)
    Vdegree=VdegreeLimit(1);
    disp(['Vertical degree should be more than ',num2str(VdegreeLimit(1))]);
end
if Vdegree>VdegreeLimit(2)
    Vdegree=VdegreeLimit(2);
    disp(['Vertical degree should be less than ',num2str(VdegreeLimit(2))]);
end

% RdegreeLimit = [-12 12];
% if Rdegree<RdegreeLimit(1)
%     Rdegree=RdegreeLimit(1);
%     disp(['Roll degree should be more than ',num2str(RdegreeLimit(1))]);
% end
% if Rdegree>RdegreeLimit(2)
%     Rdegree=RdegreeLimit(2);
%     disp(['Roll degree should be less than ',num2str(RdegreeLimit(2))]);
% end

%% Instrument Connection
% Find a serial port object.
obj1 = instrfind('Type', 'serial', 'Port', COMNum, 'Tag', '');%Change the port # here

% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = serial(COMNum);
else
    fclose(obj1);
    obj1 = obj1(1);
end
% Configure instrument object, obj1.
set(obj1, 'Timeout', 0.1);
set(obj1, 'Terminator', {'CR/LF','CR/LF'});
set(obj1, 'BaudRate', 57600);

% Connect to instrument object, obj1.
fopen(obj1);

%% Instrument Configuration and Control
% Communicating with instrument object, obj1.
fprintf(obj1, 'GG*'); %Initialization
% fprintf(obj1, 'GGN+0cz00');% Turn off home locatoin, enable confirmation response
fprintf(obj1, 'X0N-0cz00');
fprintf(obj1, 'Y0N-0cz00');
fprintf(obj1, 'Z0N-0cz00');
% fprintf(obj1, 'X0F');%Full Step mode
%fgets(obj1)
%% Step Mode
fprintf(obj1, 'X0H4');%Half Step mode
fprintf(obj1, 'Y0H4');%Half Step mode
fprintf(obj1, 'Z0H4');%Half Step mode
% X0H1 – Full Stepping 
% X0H2 – 1/2 Stepping 
% X0H3 – 1/8 microstepping 
% X0H4 – 1/16 microstepping 
% X0H5 – 1/32 microstepping
% fscanf(obj1)
%% Amps Active current
% fprintf(obj1, 'X0P64,15,0');
fprintf(obj1, 'X0P3,200,75,10');
fprintf(obj1, 'Y0P3,200,75,10');
fprintf(obj1, 'Z0P3,200,75,10');
% fprintf(obj1, 'X0P64,0,0');
% X0P64,15,0 – .37 Amps Active current, .09 Amps holding current 
% X0P64,0,0 – .37 Amps Active current,  NO holding current 


fprintf(obj1, 'X0B1600');%Begin velocity (10-300)
fprintf(obj1, 'X0E1600');% End velocity (20-2000)
fprintf(obj1, 'Y0B1600');%Begin velocity (10-300)
fprintf(obj1, 'Y0E1600');% End velocity (20-2000)
fprintf(obj1, 'Z0B1600');%Begin velocity (10-300)
fprintf(obj1, 'Z0E1600');% End velocity (20-2000)

fprintf(obj1, 'X0S1');% Acceleration slope (1-200)
fprintf(obj1, 'Y0S1');% Acceleration slope (1-200)
fprintf(obj1, 'Z0S1');% Acceleration slope (1-200)
% fscanf(obj1)
Hdeg_per_step = 16*8*5;% steps per degreeS
% degree = 20;
Hsteps = Hdeg_per_step * Hdegree;

if Hsteps>0
    fprintf(obj1, ['X0RYY+',num2str(Hsteps)]);
else
    fprintf(obj1, ['X0RYY-',num2str(-Hsteps)]);
end

Vdeg_per_step = 16*8*5;% steps per degreeS
% degree = 20;
Vsteps = Vdeg_per_step * Vdegree;
if Vsteps>0
    fprintf(obj1, ['Y0RYY+',num2str(Vsteps)]);
else
    fprintf(obj1, ['Y0RYY-',num2str(-Vsteps)]);
end


Rdeg_per_step = 16*10*5;% steps per degreeS
% degree = 20;
Rsteps = Rdeg_per_step * Rdegree;
if Rsteps>0
    fprintf(obj1, ['Z0RYY+',num2str(Rsteps)]);
else
    fprintf(obj1, ['Z0RYY-',num2str(-Rsteps)]);
end

%% Get current position
if(0)
    fprintf(obj1, 'X0m');
    position = fscanf(obj1);
    position = str2double(position(4:end))/deg_per_step;
    
    fprintf(obj1, 'Y0m');
    position = fscanf(obj1);
    position = str2double(position(4:end))/deg_per_step;
end
% fprintf(obj1, 'X0M(0)');

% fprintf(obj1, 'X0m');
% X0M(x)- Move to absolute position (Steps +/-) 
% X0I1- Read Home switch Status (1 not tripped , 0 Tripped) 
% X0I3- Read Max switch  Status (1 not tripped , 0 Tripped) 
% X0m- Read current position 
%% Wait until rotation finished
% if(obj1.BytesAvailable==0)
%    error('Error: Rotator is not working!Check Connection!'); 
% end
n = 0;
while(obj1.BytesAvailable~=96)
%     fscanf(obj1)
    if rem(n,10000) ==0
        fprintf('Rotating...\n');
    end
    n = n + 1;
end
%% Disconnect and Clean Up
% Disconnect from instrument object, obj1.
fprintf('Rotation Complete!\n');
% fscanf(obj1)
%% Get current position
if(0)
    fprintf(obj1, 'X0m');
    position = fscanf(obj1);
    position = str2double(position(4:end))/Hdeg_per_step;
    
    fprintf(obj1, 'Y0m');
    position = fscanf(obj1);
    position = str2double(position(4:end))/Vdeg_per_step;
end
%%
fclose(obj1);

% end
%{
Command List Replace X with Y for other axis  
 GG* - Initialization Command – MUST BE ISSUED  
X0N-0cz00 –Initialization Command- MUST BE ISSUED , sets both home/limit switches to normally open 
X0Hx – Step Mode – See section in manual for options X0P64,15,0 – .37 Amps Active current, .09 Amps holding current X0P64,0,0 – .37 Amps Active current,  NO holding current X0B(x) - Begin Velocity Range 10-30,000 (steps/s)   X0E(x) - End Velocity Range 20-30,000  (steps/s)   X0S(x) - Acceleration slope Range 1-200 (1-3 recommended)  X0RYY(x) (x) - Relative Steps Command + or – followed by number of steps, monitor home/limit sensors X0M(x)- Move to absolute position (Steps +/-) X0I1- Read Home switch Status (1 not tripped , 0 Tripped) X0I3- Read Max switch  Status (1 not tripped , 0 Tripped) X0m- Read current position 
%}
