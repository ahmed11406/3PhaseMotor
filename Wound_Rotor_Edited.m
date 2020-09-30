% Author: Ahmed Adel Hassan Elmashhor, Section: 1.
% Subject: Plotting the torque-speed curves of the wound rotor IM when
% external resistance is added to its rotor.
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
R2org= input('Rotor resistance per phase in Ohms: \n ');
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

 % initialization of developed torque vector that contains the points to be plotted.
 
 Tdev=zeros(1,1000); 
 
 % Step 4) Calculating the developed torque  for each case using nested loops:
fprintf('Enter the variable-resistor nominal step in per-unit Ohms/step:\n(Ex. enter the value "0.1" p.u Ohms/step given rated rotor resistance at 1 p.u)');
disp('Please note that the program draws a plot per each step in the range');
Step=input('');
disp('Enter the maximum range of the external variable-resistor in p.u rated:');
disp('For example, enter "5" for a max. range equal 5 times the rated rotor resistance.');
Range=input('');
 for a = 0:Step:Range 
     R2= R2org + a * R2org; 
    for s = 1:1:1000
Zeq = complex(R1,X1) + complex((-Xm*X2*Rc),(Xm*(R2/(s/1000))*Rc))/(complex(((R2/(s/1000))*Rc-X2*Xm),(Rc*X2+ Xm*((R2/(s/1000))+Rc))));
I1 = Vph/Zeq;
I2 = I1 * complex(0,Xm*Rc)/(complex(((R2/(s/1000))*Rc-X2*Xm),(Rc*X2+ Xm*((R2/(s/1000))+Rc))));
Pgap = 3 * (abs(I2)^2)*(R2/(s/1000));
Tdev(1001-s) = Pgap/ws;
    end
s = 1:-0.001:0.001;
Marker = ['Rex = ',num2str(a),' R2']; % Variable-text vector that contains the curves annotations.
plot(s,Tdev,'DisplayName',Marker)     % Assign each curve to its annotation.
hold on                               % wait for another curve to be added to the same plot.
grid on                               % make the plot screen gridded. 
set(gca, 'XDir','reverse')            % Reverses the direction of X-axis (the slip drown from 1 to 0 instead of 0 to 1).
legend show                           % show the legend on the plot window.
end
return;                               % exits the script.