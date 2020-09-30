% Author: Ahmed Adel Hassan Elmashhor, Section: 1.
% Subject: Plotting the torque-speed curves of the squirrel-cage IM in the
% following three cases:
% A) Variable voltage, constant frequency supply
% B) constant voltage, variable frequency supply
% C) Variable voltage, variable frequency supply (V/F control)
%Note: in this script, the author uses the exact equivalent circuit model of the
%three phase induction motor for all calculations.

clear clcs; % clears all current workplace variables
% Step 1) Getting the name-plate data of the motor entered by the user:

fprintf('Please! Enter the following information about the machine reffered to stator:\n \n')
Fs= input('The rated frequency in Hz: \n '); 
V= input('The rated stator line voltage in Volts: \n ');
Connection= input('Select the connection type:\nEnter 1 for star type, or\nEnter 2 for delta type.\n ');
P= input('The number of poles: \n ');  
R1= input('Stator resistance per phase in Ohms: \n '); 
R2= input('Rotor resistance per phase in Ohms: \n ');
X1= input('Stator reactance per phase in Ohms: \n ');
X2= input('Rotor reactance per phase in Ohms: \n ');
Xm= input('Magnetization reactance per phase in Ohms: \n ');
Rc= input('Core resistance per phase in Ohms: \n ');

% Step 2) Calculating the W_sync and ns:

ns = 120*Fs/P;
ws = 2*pi*ns/60;

% Step 3) Calculating V per phase:

if (Connection==1) % Star.
        Vph = V/sqrt(3);
end
if (Connection==2) % Delta.
        Vph = V; 
end

% Step 4) Enabling the user to select certain case to plot:

 fprintf('\nSelect which case to plot:\n');   
 disp('Enter 1 for Case A)Variable voltage, constant frequency supply.');
 disp('Enter 2 for Case B)Constant voltage, variable frequency supply.');
 disp('Enter 3 for Case C)Variable voltage, variable frequency supply (V/F control).');
 Case_Selector=input('');
 
 % initialization of developed torque vector that contains the points to be plotted.
 
 Tdev=zeros(1,1000); 
 nr=zeros(1,1000);
 
 % Step 5) Calculating the developed torque  for each case using nested loops: 
 
 if Case_Selector==1 
% Case A) Variable voltage, constant frequency supply.
%--------------------------------------------------------------------------
fprintf('Enter the voltage controller nominal step in per-unit volt/step:\n(Ex. enter the value "0.1" p.u volt/step given rated voltage at 1 p.u)');
disp('Please note that the program draws a plot per each voltage step');
Step=input('');
% Let the supply voltage varies over a range starting from 0 to the rated value with a user percentage step.
for a = 0:Step:1 
    for s = 1:1:1000
Zeq = complex(R1,X1) + complex((-Xm*X2*Rc),(Xm*(R2/(s/1000))*Rc))/(complex(((R2/(s/1000))*Rc-X2*Xm),(Rc*X2+ Xm*((R2/(s/1000))+Rc))));
I1 = a * Vph/Zeq;
I2 = I1 * complex(0,Xm*Rc)/(complex(((R2/(s/1000))*Rc-X2*Xm),(Rc*X2+ Xm*((R2/(s/1000))+Rc))));
Pgap = 3 * (abs(I2)^2)*(R2/(s/1000));
Tdev(1001-s) = Pgap/ws;
    end
s = 1:-0.001:0.001;
Marker = ['V= ',num2str(a),' V_rated'];    % Variable-text vector that contains the curves annotations.
plot(s,Tdev,'DisplayName',Marker)          % Assign each curve to its annotation.
hold on
grid on
set(gca, 'XDir','reverse')
legend show
end
%--------------------------------------------------------------------------
 end
 
 if Case_Selector==2
% Case B) constant voltage, variable frequency supply.
%--------------------------------------------------------------------------
fprintf('Enter the frequency controller nominal step in per-unit Hz/step:\n(Ex. enter the value "0.1" p.u Hz/step given rated frequency at 1 p.u)\n');
disp('Please note that the program draws a plot per each frequency step');
Step=input('');
disp('Enter the maximum range of the frequency in p.u:');
disp('For example, enter "5" for a max. range equal 5 times the rated frequency.');
Range=input('');
% Let the supply frequency varies over a range starting from 0 to Range times the rated value with a percentage step of Step.
for a = 0:Step:Range
    for s = 1:1:1000
Zeq = complex(R1,a*X1) + complex((-a*Xm*a*X2*Rc),(a*Xm*(R2/(s/1000))*Rc))/(complex(((R2/(s/1000))*Rc-a*X2*a*Xm),(Rc*a*X2+ a*Xm*((R2/(s/1000))+Rc))));
I1 = Vph/Zeq;
I2 = I1 * complex(0,a*Xm*Rc)/(complex(((R2/(s/1000))*Rc-a*X2*a*Xm),(Rc*a*X2+ a*Xm*((R2/(s/1000))+Rc))));
Pgap = 3 *(abs(I2)^2)*(R2/(s/1000));
Tdev(1001-s) = Pgap/a*ws;
ns = a*120*Fs/P;
nr(1001-s)= (1-(s/1000))*ns;
    end
Marker = ['F= ',num2str(a),' F-rated'];    % Variable-text vector that contains the curves annotations.
plot(nr,Tdev,'DisplayName',Marker)         % Assign each curve to its annotation.
hold on
grid on
legend show
end
%--------------------------------------------------------------------------
 end

 if Case_Selector==3
% Case C) Variable voltage, variable frequency supply (V/F control).
%--------------------------------------------------------------------------
fprintf('Enter the V/F controller nominal step in per-unit Volt/Hz/step:\n(Ex. enter the value "0.1" p.u Volt/Hz/step given rated voltage at 1 p.u)\n');
disp('Please note that the program draws a plot per each step');
Step=input('');
% Let the supply frequency and voltage varies over a range starting from 0 to the voltage rated value with a percentage step of Step.
for a = 0:Step:1
    for s = 1:1:1000
Zeq = complex(R1,a*X1) + complex((-a*Xm*a*X2*Rc),(a*Xm*(R2/(s/1000))*Rc))/(complex(((R2/(s/1000))*Rc-a*X2*a*Xm),(Rc*a*X2+ a*Xm*((R2/(s/1000))+Rc))));
I1 = a*Vph/Zeq;
I2 = I1 * complex(0,a*Xm*Rc)/(complex(((R2/(s/1000))*Rc-a*X2*a*Xm),(Rc*a*X2+ a*Xm*((R2/(s/1000))+Rc))));
Pgap = 3 *(abs(I2)^2)*(R2/(s/1000));
Tdev(1001-s) = Pgap/a*ws;
ns = a*120*Fs/P;
nr(1001-s)= (1-(s/1000))*ns;
    end
Marker = ['V/F ratio ',num2str(a),' p.u']; % Variable-text vector that contains the curves annotations.
plot(nr,Tdev,'DisplayName',Marker)         % Assign each curve to its annotation.
hold on
grid on
legend show
end
%--------------------------------------------------------------------------
 end
return;