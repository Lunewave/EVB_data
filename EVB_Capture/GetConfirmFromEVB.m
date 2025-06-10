% Get confirmation String from EVB
function Confirmed = GetConfirmFromEVB(evb, string)

readback = ['init'];
pause(1);

bytesAvailable= evb.NumBytesAvailable;
if (bytesAvailable~=0)
    readback= read(evb, bytesAvailable, "string");
end
read_index=1;

while (contains(readback, string) ~= 1)
    pause(0.1);

    bytesAvailable= evb.NumBytesAvailable;
    if (bytesAvailable~=0)
        readback= read(evb, bytesAvailable, "string");
    end
   
    read_index=read_index+1;
    if read_index > 200
        disp('readback error')
        beep; pause(1);
        beep; pause(1);
        beep; pause(1);
        pause;
        continue;
    end
end

Confirmed = 1;

end
