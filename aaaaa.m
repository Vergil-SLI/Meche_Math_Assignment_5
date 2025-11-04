function aaaaa
x = 1;
y = 3;
theta = 0;
Plist_box = [1,1,1;1,1,1];

box_params = struct();

box_params.m =  ;%mass of the box
box_params.I =  %moment of inertia w/respect to centroid
box_params.g =  %acceleration due to gravity
box_params.k_list = %list of spring stiffnesses
box_params.l0_list = %list of spring natural lengths

Plist_world = compute_rbt(x,y,theta,Plist_box)

box_params.P_world =  Plist_world; %2 x n list of static mounting points for the spring (in the world frame)
box_params.P_box =  Plist_box; %2 x n list of mounting pointsfor the spring (in the box frame)


[ax,ay,atheta] = compute_accel(x,y,theta,box_params)

end