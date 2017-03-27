function [m_out,r_out, I_out, rotor_locations,J_out, cds_out, A_cross_out,rotor_dist_out,wdir_out,cls_out, cpqr_out] = get_dimensions

global I_r;
global r;
global r_loc;
global m;
global J;
global cds;
global A_cross;
global rotor_dist;
global w_dir;
global cls;
global c_pqr;
m_out = m;
I_out = I_r;
r_out = r;
rotor_locations = r_loc;
J_out = J;
cds_out = cds;
A_cross_out = A_cross;
rotor_dist_out=rotor_dist;
wdir_out=w_dir;
cls_out=cls;
cpqr_out=c_pqr;

end

