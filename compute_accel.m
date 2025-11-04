%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of mounting point Qs (in the world frame)
%box_params.P_box: 2 x n list of mounting points Ps (in the box frame)
%
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)

    F_sum = 0;
    Q = compute_rbt(x,y,theta,box_params.P_box); %P_world = Q

    for i = 1:length(box_params.k_list)
        PA = box_params.P_box(1,i);
        PB = box_params.P_box(2,i);
        F = compute_spring_force(box_params.k_list(i),box_params.l0_list(i),PA,PB);
        F_sum = F_sum + F;
    end

end