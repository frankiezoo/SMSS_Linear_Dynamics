syms rho_fe_x rho_fe_y rho_fe_z real
syms rho_me_x rho_me_y rho_me_z real
syms dr_x dr_y dr_z real
syms m_x m_y m_z dm_x dm_y dm_z ds real
syms m_fx m_fy m_fz m_mx m_my m_mz real
syms mu_0 mu a_sx a_sy a_sz real

% superconductor normal
a_s=[a_sx a_sy a_sz]';
% equilibrium position
rho_fe = [rho_fe_x rho_fe_y rho_fe_z]';
rho_me = [rho_me_x rho_me_y rho_me_z]';
dr = [dr_x dr_y dr_z]';
drm = 2*(a_s*a_s')*dr;
% equilibrium magnet orientation
m_e = [m_x m_y m_z]';
dm = [dm_x dm_y dm_z]';
dq=[dm;ds];
% equilibrium image orientation
m_fe = [m_fx m_fy m_fz]';
m_me = [m_mx m_my m_mz]';

var = [dr_x dr_y dr_z dm_x dm_y dm_z ds]';

% cross terms for convenience
rho_fex=crs(rho_fe);
rho_mex=crs(rho_me);
m_ex=crs(m_e);
m_fex=crs(m_fe);
m_mex=crs(m_me);

%% derive matrix terms for force
dF_f_lin_drnorm=3*mu_0/(4*pi)*...
    -5/norm(rho_fe)^7*(rho_fe*m_e'*m_fe + m_e*rho_fe'*m_fe + m_fe*rho_fe'*m_e ...
        -5/norm(rho_fe)^2*(rho_fe*rho_fe')*m_e*rho_fe'*m_fe)*rho_fe'*dr;
    
dF_f_lin_drmat=3*mu_0/(4*pi*norm(rho_fe)^5)*(...
    (m_fex*m_ex+m_ex*m_fex-2*m_e'*m_fe*eye(3)-5/norm(rho_fe)^2*(rho_fe*(rho_fex*m_fe)'*m_ex...
    -rho_fe*(rho_fex*m_e)'*m_fex+(rho_fex*m_e)'*(rho_fex*m_fe)))*dr+...
    (-m_fex*rho_fex+crs(rho_fex*m_fe)-2*rho_fe*m_fe'+5/norm(rho_fe)^2*rho_fe*(rho_fex*m_fe)'*rho_fex)*[mu*eye(3) zeros(3,1)]*dq);


[F_f b_f]=equationsToMatrix(dF_f_lin_drnorm+dF_f_lin_drmat,var)

dF_m_lin_drnorm=3*mu_0/(4*pi)*...
    -5/norm(rho_me)^7*(rho_me*m_e'*m_me + m_e*rho_me'*m_me + m_me*rho_me'*m_e ...
        -5/norm(rho_me)^2*(rho_me*rho_me')*m_e*rho_me'*m_me)*rho_me'*2*(a_s*a_s')*dr;
    
dF_m_lin_drmat=3*mu_0/(4*pi*norm(rho_me)^5)*(...
    (m_mex*m_ex+m_ex*m_mex-2*m_e'*m_me*eye(3)-5/norm(rho_me)^2*(rho_me*(rho_mex*m_me)'*m_ex...
    -rho_me*(rho_mex*m_e)'*m_mex+(rho_mex*m_e)'*(rho_mex*m_me)))*(2*(a_s*a_s'))*dr+...
    (-m_mex*rho_mex+crs(rho_mex*m_me)-2*rho_me*m_me'+5/norm(rho_me)^2*rho_me*(rho_mex*m_me)'*rho_mex...
    +(eye(3)-2*(a_s*a_s'))*(crs(rho_mex*m_e)-m_ex*rho_mex-2*rho_me*m_e'+5/norm(rho_me)^2*rho_me*(rho_mex*m_e)'*rho_mex))*[mu*eye(3) zeros(3,1)]*dq);
  
[F_m b_m]=equationsToMatrix(dF_m_lin_drnorm+dF_m_lin_drmat,var)

A_force=F_f+F_m;

%% derive matrix terms for torque
dtau_f_lin_drnorm=mu_0/(4*pi)*...
    -3/norm(rho_fe)^5*(3/norm(rho_fe)^2*(m_e'*rho_fe)*(m_fex*rho_fe)+m_ex*m_fe)*rho_fe'*dr;
    
dtau_f_lin_drmat=mu_0/(4*pi*norm(rho_fe)^3)*(...
    3/norm(rho_fe)^2*(m_e'*rho_fe*m_fex+m_fex*rho_fe*m_e')*dr+...
    (3/norm(rho_fe)^2*(m_fex*(rho_fe*rho_fe'))-m_fex)*[mu*eye(3) zeros(3,1)]*dq);

[tau_f b_f]=equationsToMatrix(dtau_f_lin_drnorm+dtau_f_lin_drmat,var)

dtau_m_lin_drnorm=mu_0/(4*pi)*...
    -3/norm(rho_me)^5*(3/norm(rho_me)^2*(m_e'*rho_me)*(m_mex*rho_me)+m_ex*m_me)*rho_me'*2*(a_s*a_s')*dr;
    
dtau_m_lin_drmat=mu_0/(4*pi*norm(rho_me)^3)*(...
    3/norm(rho_me)^2*(m_e'*rho_me*m_mex+m_mex*rho_me*m_e')*2*(a_s*a_s')*dr+...
    (3/norm(rho_me)^2*(-m_e'*rho_me*rho_mex)+m_ex)*(eye(3)-2*(a_s*a_s'))*[mu*eye(3) zeros(3,1)]*dq+...
    (3/norm(rho_me)^2*(m_mex*(rho_me*rho_me'))-m_mex)*[mu*eye(3) zeros(3,1)]*dq);

[tau_m b_m]=equationsToMatrix(dtau_m_lin_drnorm+dtau_m_lin_drmat,var)

A_torque=tau_f+tau_m;

%% Complete Jacobian [dF/dr dF/dm; dtau/dr dtau/dm]
dFtau_drm=simplify([A_force;A_torque])