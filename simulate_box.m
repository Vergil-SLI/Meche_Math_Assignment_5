function simulate_box()
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
    rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    x0 = 0; 
    y0 = 0;
    theta0 = pi/6;
    vx0 = 0; 
    vy0 = 1; 
    vtheta0 = pi/4; 

    V0 = [x0; y0; theta0; vx0; vy0; vtheta0];
    tspan = [0; 10];

    % define values of RK method (explicit midpoint)
    h_ref = 0.01;
    BT_struct = struct();
    BT_struct.A = [0, 0; 0.5, 0]; % matrix of a_{ij} values
    BT_struct.B = [0, 1];% vector of b_i values
    BT_struct.C = [0, 0.5]; % vector of c_i values
    
    % run the integration
    % Vlist is a list of [x;y;theta;dxdt;dydt;dthetadt]
    [tlist,Vlist] = explicit_RK_fixed_step_integration(rate_func,tspan,V0,h_ref,BT_struct);

    % plot the springs as animation
    num_zigs = 5;
    w = .1;
    hold on;
    spring_plot_1 = initialize_spring_plot(num_zigs,w);
    spring_plot_2 = initialize_spring_plot(num_zigs,w);
    spring_plot_3 = initialize_spring_plot(num_zigs,w);
    spring_plot_4 = initialize_spring_plot(num_zigs,w);
    box_plot = plot(0, 0, "black");
    axis equal; axis square;
    axis([-10,10,-10,10])


    for i = 1:length(tlist)
        % find the springs' endpoint that's attached to box
        x = Vlist(1, i);
        y = Vlist(2, i);
        theta = Vlist(3, i);
        P_world = compute_rbt(x,y,theta,P_box);

        % draw springs
        update_spring_plot(spring_plot_1,P_world(:, 1),Q_world(:, 1))
        update_spring_plot(spring_plot_2,P_world(:, 2),Q_world(:, 2))
        update_spring_plot(spring_plot_3,P_world(:, 3),Q_world(:, 3))
        update_spring_plot(spring_plot_4,P_world(:, 4),Q_world(:, 4))

        % find box coord & plot it
        box_boundary_world = compute_rbt(x, y, theta, boundary_pts);
        set(box_plot, 'xdata', box_boundary_world(1, :), 'ydata', box_boundary_world(2, :));
        drawnow;
    end
    
end


%updates spring plotting object so that spring is plotted
%with ends located at points P1 and P2
function update_spring_plot(spring_plot_struct,P1,P2)
    dP = P2-P1;
    R = [dP(1),-dP(2)/norm(dP);dP(2),dP(1)/norm(dP)];
    plot_pts = R*spring_plot_struct.zig_zag;
    set(spring_plot_struct.line_plot,...
    'xdata',plot_pts(1,:)+P1(1),...
    'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,...
    'xdata',[P1(1),P2(1)],...
    'ydata',[P1(2),P2(2)]);
end


%create a struct containing plotting info for a single spring
%INPUTS:
%num_zigs: number of zig zags in spring drawing
%w: width of the spring drawing
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
    spring_plot_struct = struct();
    zig_ending = [.25,.75,1; ...
    -1,1,0];
    zig_zag = zeros(2,3+3*num_zigs);
    zig_zag(:,1) = [-.5;0];
    zig_zag(:,end) = [num_zigs+.5;0];
    
    for n = 0:(num_zigs-1)
        zig_zag(:,(3+3*n):2+3*(n+1)) = zig_ending + [n,n,n;0,0,0];
    end
    
    zig_zag(1,:)=(zig_zag(1,:)-zig_zag(1,1))/(zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:)=zig_zag(2,:)*w;
    spring_plot_struct.zig_zag = zig_zag;
    spring_plot_struct.line_plot = plot(0,0,'k','linewidth',2);
    spring_plot_struct.point_plot = plot(0,0,'ro','markerfacecolor','r','markersize',7);
end