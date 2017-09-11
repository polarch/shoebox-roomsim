function h_ax = plot_scene(room, src, rec, h_ax)
%PLOT_SCENE Simple plot of sources-receivers-room geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT_SCENE.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4, h_ax = axes;
else axes(h_ax), end

% room dimensions
l = room(1);
w = room(2);
h = room(3);

% source coordinates
S = src;
% move origin to the centre of the room
S_x = S(:, 1) - l/2;
S_y = w/2 - S(:, 2);
S_z = S(:, 3) - h/2;

% receiver coordinates
R = rec;
% move origin to the centre of the room
R_x = R(:, 1) - l/2;
R_y = w/2 - R(:, 2);
R_z = R(: ,3) - h/2;

room = makeRectangle(l, w, h);

camlight
box on
view(30,30)
set(gca,'xlim',[-1.5*l 1.5*l],'ylim',[-1.5*w 1.5*w],...
    'zlim',[-1.5*h 1.5*h])
hold on

patchFace(room, 'facecolor','g', 'FaceAlpha',0.1);

%plot 3d axes
line([0 l],[0 0], [0 0],'color','r');
text(l,0,0,'x','Color','r','FontSize',24);
line([0 0],[0 w], [0 0],'color','g');
text(0,w,0,'y','Color','g','FontSize',24);
line([0 0],[0 0], [0 h],'color','b');
text(0,0,h,'z','Color','b','FontSize',24);

% set up unit sphere information
numSphereFaces = 10;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);

spheresRadius = 0.2;

for i=1:size(S, 1)
    sphereX = S_x(i) + unitSphereX*spheresRadius;
    sphereY = S_y(i) + unitSphereY*spheresRadius;
    sphereZ = S_z(i) + unitSphereZ*spheresRadius;

    surface(sphereX, sphereY, sphereZ,'FaceColor', 'b');            
end

for i=1:size(R, 1)
    sphereX = R_x(i) + unitSphereX*spheresRadius;
    sphereY = R_y(i) + unitSphereY*spheresRadius;
    sphereZ = R_z(i) + unitSphereZ*spheresRadius;

    surface(sphereX, sphereY, sphereZ,'FaceColor', 'r');            
end

light('Position',[2*l 0 2*h])        
axis equal

end

function box = makeRectangle(l, w, h)

% get vertices, starting form down left corner, length l along x
% dimension, width w along y dimension, height along z
%
%      ^ x      1.---.4
%      |         |   |
%      |         |   |
% y <--.        2.---.3
%

vertices(1,:) = [l/2, w/2, -h/2];
vertices(2,:) = [-l/2, w/2, -h/2];
vertices(3,:) = [-l/2, -w/2, -h/2];
vertices(4,:) = [l/2, -w/2, -h/2];

vertices(5,:) = [l/2, w/2, h/2];
vertices(6,:) = [-l/2, w/2, h/2];
vertices(7,:) = [-l/2, -w/2, h/2];
vertices(8,:) = [l/2, -w/2, h/2];

box.vert = vertices;
box.face(1,:) = [4 3 2 1];
box.face(2,:) = [5 6 7 8];
box.face(3,:) = [5 8 4 1];
box.face(4,:) = [5 6 2 1];
box.face(5,:) = [6 7 3 2];
box.face(6,:) = [7 8 4 3];

end

function patchFace(poly,varargin)
%PATCHFACE Plots a sequence of polygon faces defined by a structure.
    if iscell(poly.face)
        for n = 1:length(poly.face)
            patch('faces', poly.face{n}, 'vertices', poly.vert, varargin{:});
        end
    elseif ismatrix(poly.face)
        for n = 1:size(poly.face, 1)
            patch('faces', poly.face(n, :), 'vertices', poly.vert, varargin{:});
        end
    else
        error('poly.face should be either a cell array or a matrix')
    end
end
