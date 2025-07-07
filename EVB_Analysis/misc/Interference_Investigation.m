D1 = 4.5; %Line of sight distance from horn antenna to Dipole in meters
h = 0.07; %height of dipole above the ground plane in meters
theta = 30; %Elevation angle down from horizontal in degrees
phi = 90-theta; %Angle of reflection off ground plane
f = 2.456; %frequency in GHz




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = deg2rad(theta);
phi = deg2rad(phi);
d2 = h/sin(theta);
c = 3*10^8;
lambda = c/(f*10^9);
a = asin(d2*sin(2*phi)/D1); %Angle between line of sight signal and reflected signal
b = pi-2*phi-a;
d1 = D1*sin(b)/sin(2*phi);

D2 = d1+d2;

diff = D2-D1;

cycle_diff = diff/lambda;

disp(['The two signal are ' num2str(cycle_diff) ' apart.'])