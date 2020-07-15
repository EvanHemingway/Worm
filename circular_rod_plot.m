function [V,rod]=circular_rod_plot(r,d1,d2,N)
    %This function takes in centerline position vectors with components x,y
    %and z at every discrete reference arc length. Additionally, it takes
    %in the directors scaled such that they point at the boundary of the
    %rod, also as functions of discrete reference arc length. It plots the
    %current configuration of the rod assuming a circular cross-section in
    %the reference configuration with the centerline chosen as the line of
    %geometric centers.
    
    %r is a CX3 array, where C is the number of discrete reference arc
    %lengths
    
    %Similarly, d1 and d2 are CX3
    
    %N is the number of points desired in each cross-section

D = size(r);
C = D(1); %Number of cross-sections

thet=linspace(0,2*pi*(1-1/N),N);
thet = thet';

V=zeros(C*N,3);
for i=1:C
    c=cos(thet)*d1(i,:)+sin(thet)*d2(i,:);
    rb=r(i,:)+c;
    V(N*(i-1)+1:N*i,:)=rb;
end

F = zeros(N*(C-1),5); %quad elements
for i=1:N*(C-1)
    if mod(i,N)==0
        F(i,:)=[i,i+1-N,i+1,i+N,i];
    else
        F(i,:)=[i,i+1,i+1+N,i+N,i];
    end
end

rod=patch('Faces',F,'Vertices',V,'FaceColor','red','EdgeColor','none');

%patch( [-10 -10 10 10] , [50 0 0 50], [-1 -1 -1 -1], [1 1 1 1])

%view(53,41)
%camzoom(0.3)
% xlim([0 1.2])
% ylim([-1 1.5])
% zlim([0 1])

% xticks([-10 0 10])
% yticks([0 10 20 30 40 50])
% zticks([0 10 20])

grid on

ax = gca;
ax.Clipping = 'off';
end