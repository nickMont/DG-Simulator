function R = calculate_R_from_euler(e)
%RBI (rotation to body from inertial)
phi = e(1);
theta = e(2);
psi = e(3);

R =[cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta)
    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta)
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];

end

