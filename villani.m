function [F,tau]=villani(mu_0,mu_mj,mu_i,rho_mji)
%all equations taken from Yung/Villani Force and Torque Between Magnetic Dipoles
%test script named Villani_test, passed all test cases in Villani force paper
mu_mj_norm=norm(mu_mj);
mu_mj_hat=mu_mj/mu_mj_norm;

mu_i_norm=norm(mu_i);
mu_i_hat=mu_i/mu_i_norm;

rho_mji_norm=norm(rho_mji);
rho_mji_hat=rho_mji/rho_mji_norm;

F=3*mu_0*mu_mj_norm*mu_i_norm/4/pi/rho_mji_norm^4*[...
    rho_mji_hat*dot(mu_mj_hat,mu_i_hat),...
    mu_mj_hat*dot(rho_mji_hat,mu_i_hat),...
    mu_i_hat*dot(rho_mji_hat,mu_mj_hat),...
    -5*rho_mji_hat*dot(rho_mji_hat,mu_mj_hat)*dot(rho_mji_hat,mu_i_hat)];
tau=mu_0*mu_mj_norm*mu_i_norm/4/pi/rho_mji_norm^3*...
    (3*dot(mu_mj_hat,rho_mji_hat)*cross(mu_i_hat,rho_mji_hat)+...
    cross(mu_mj_hat,mu_i_hat));

end