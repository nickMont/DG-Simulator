%Intialization of quad(s)

numquads=1;

e_w_x_v_prev(:,1,1) = [0;0;0; 0;0;0; 0;0;0 ; 0;0;0];
e_w_x_v_des = zeros(12,numquads);
e_w_x_v_des(:,1)=[0;0;0; 0;0;0; 10;0;0; 0;0;0];

% waypoints=[e_w_x_v_des1 e_w_x_v_des2];  %
% %set_destination_point(waypoints(:,1));


pos_close=false;
angle_close=false;
waypoint_number=1;
%[~,max_waypoints]=size(waypoints);
%e_w_x_v_des = waypoints(:,waypoint_number);

F_control_max=50;
F_control_min=-10;
min_thrust=2;
max_thrust=60;
set_f_bounds([F_control_min;F_control_max],[min_thrust;max_thrust]);

%Moments of inertia about various axes
Ixx = .75;
Iyy = .8;
Izz = 1;
Ixy = 0;
Ixz = 0;
Iyz = 0;
J=ones(3,3,numquads);
J(:,:,1) = [Ixx Ixy Ixz
    Ixy Iyy Iyz
    Ixz Iyz Izz];
if numquads>=2
    for i=2:numquads
        J(:,:,i)=J(:,:,1);
    end
end

m = 1.5*ones(1,numquads); %kg

%Physical properties of rotors (all assumed identical)
I_rotor = .02*ones(1,numquads); %moment of inertia of a single rotor
r_rotor = .25*ones(1,numquads);
r_dis = 1; %distance between cg of quad and rotor center

cds = [.1;.1;.2]*ones(1,numquads); %Drag coefficients in body frame for velocity directed along each axis
%ie cds(1) = cd_x, or the drag coefficient corresponding to a velocity in x
cls = [.01;.01]*ones(1,numquads);
c_pqr=[0;0;0]*ones(1,numquads);

cross_area = [.2;.2,;.5]*ones(1,numquads);
%Projected cross-sectional area from a viewpoint along each axis

%Rotor locations as vectors in body frame:
%    ^      (nose)
% 4     1
%    x
% 3     2
%                   |
% . is origin, -x-, y 
%                   |
r1loc = r_dis*unit_vector([1;1;0]);
r2loc = r_dis*unit_vector([1;-1;0]);
r3loc = r_dis*unit_vector([-1;-1;0]);
r4loc = r_dis*unit_vector([-1;1;0]);
rotor_loc=zeros(3,4,numquads);
rotor_loc(:,:,1)=[r1loc r2loc r3loc r4loc];
if numquads>=2
    for i=2:numquads
        rotor_loc(:,:,i)=rotor_loc(:,:,1);
    end
end

wdir=[1 -1 1 -1]'*ones(1,numquads);

set_dimensions(m, I_rotor, r_rotor, rotor_loc, J, cds, cross_area,r_dis,wdir,cls,c_pqr);
%NOTE how axes are rendered with y going from +1 to -1.  It can be
%disorienting, make sure to check coordinates of engine in the plot before
%calling something as incorrect


%Set motor properties
%~,~,~,kT
set_Motor_Properties(.7,.5,.7,1.5)
%See code to see which call corresponds to which coefficient/update later

Ts = 10; %constant of how quickly wind speed changes
Tv = .2; %constant of how quickly wind direction changes
s0 = 2; %base windspeed
d0 = [.9;.1;0]; %unit vector of base wind direction
set_Air_Properties(1.225,Ts,Tv,s0,d0)