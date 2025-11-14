function equilibrium_solver()
    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
    l0 = 1.5*LG;
    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;
    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];
    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];
    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];
    Q_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];
    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];
    
    %define system parameters
    t_in = 0; % this doesn't matter, I called it!! -Sherry
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(Q_world,2));
    box_params.Q_world = Q_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;

    %load the system parameters into the rate function via an anonymous function
    rate_func = @(V_in) box_rate_func(t_in,V_in,box_params);
    x0 = 0; 
    y0 = 0;
    theta0 = pi/6;
    vx0 = 0; 
    vy0 = 1; 
    vtheta0 = pi/4; 
    V0 = [x0; y0; theta0; vx0; vy0; vtheta0];

    solver_params = struct();
    [V_eq, exit_flag] = multi_newton(rate_func,V0,solver_params);

    %Linearization________________________________________________________
    % initial values (near equilibrium)
    V0 = V_eq;
    tspan = [0; 10];
    J_approx = approximate_jacobian_for_newton(rate_func, V_eq);
    J_approx = J_approx(4:6, 1:3);

    % define values of RK method (explicit midpoint)
    h_ref = 0.01;
    BT_struct = struct();
    BT_struct.A = [0, 0; 0.5, 0]; % matrix of a_{ij} values
    BT_struct.B = [0, 1];% vector of b_i values
    BT_struct.C = [0, 0.5]; % vector of c_i values

    % the 2 functions we're comparing
    nonlinear_func = @(t,V) box_rate_func(t,V,box_params);
    linear_func = @(t,V) J_approx*(V-V_eq);

    
end


% Numerical Approximation of the Jacobian
function J = approximate_jacobian_for_newton(fun,x)
    J = [];
    h = 1e-6;

    for j = 1:length(x)
        basis_j = zeros(length(x), 1);
        basis_j(j) = 1;
        column = (fun(x + h*basis_j) - fun(x - h*basis_j)) / (2*h);
        J = [J, column];
    end
end