%Computes the linear and angular acceleration of the box
%given its current position and orientation

%INPUTS:
%x: current x position of the box's centroid
%y: current y position of the box's centroid
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.Q_world: 2 x n list of mounting point Qs (in the world frame),
%previously referred to as P_world
%box_params.P_box: 2 x n list of mounting points Ps (in the box frame)

%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
    % RENAME STRUCT NAME FROM "P_world" TO "Q_world" IF GIVEN VALS
    F_sum = 0;
    torque_sum = 0;
    P_world = compute_rbt(x,y,theta,box_params.P_box); % P vals in world frame

    for i = 1:length(box_params.k_list)
        Qi = box_params.Q_world(:,i);
        Pi = P_world(:,i);

        F = compute_spring_force(box_params.k_list(i),box_params.l0_list(i),Qi,Pi);
        torque = norm(cross([Pi; 0] - [x; y; 0], [F; 0]));

        F_sum = F_sum + F;
        torque_sum = torque_sum + torque;
    end

    trans_accel = ([0; -1*box_params.m*box_params.g] + F_sum) ./ box_params.m;
    ax = trans_accel(1);
    ay = trans_accel(2);
    atheta = torque_sum / box_params.I;
end