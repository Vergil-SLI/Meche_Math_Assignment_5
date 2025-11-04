%Computes the rigid body transformation that maps a set
%of points in the box-frame to their corresponding
%world-frame coordinates
%INPUTS:
%x: the x position of the centroid of the box
%y: the y position of the centroid of the box
%theta: the orientation of the box
%Plist_box: a 2 x n matrix of points in the box frame
%OUTPUTS:
%Plist_world: a 2 x n matrix of points describing
%the world-frame coordinates of the points in Plist_box
function Plist_world = compute_rbt(x,y,theta,Plist_box)
    Centroid = ones(2, length(Plist_box(1, :)));
    Centroid(1,:) = Centroid(1,:) * x;
    Centroid(2,:) = Centroid(2,:) * y;
    Rotation = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    
    Plist_world =(Rotation*Plist_box)+Centroid;
end