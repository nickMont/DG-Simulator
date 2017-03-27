function set_Motor_Properties(eta_in, f_in, Kq_in,Kt_in)

%eta, f, Kt, and Ka are constants
%eta--motor efficiency

%figure of merit of propeller, estimated value of 0.50 which is
%reasonable based on this paper:
%http://m-selig.ae.illinois.edu/pubs/DetersSelig-2008-AIAA-2008-6246-MicroProps.pdf

%Kq--motor coefficient Kq, see:
%http://ieeexplore.ieee.org.ezproxy.lib.utexas.edu/stamp/stamp.jsp?tp=&arnumber=4505621&tag=1

%Kt--constant of proportionality of propeller torque and propeller thrust
%as Torque = Kt * Thrust

global eta;
global figure_of_merit;
global Kq;
global Kt;

eta = eta_in;
figure_of_merit = f_in;
Kq = Kq_in;
Kt = Kt_in;


end

