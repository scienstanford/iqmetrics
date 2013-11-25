function plot_quadtree(QT,Theta, M, plot_type)

% plot_quadtree - graphical display of a quatree
%
%   plot_quadtree(QT,Theta, M, plot_type);
%
%   M is an arbitrary background image
%   QT and Theta are computed using the function 
%        [QT,Theta] = compute_quadtree(M,T,j_min,j_max,s);
%
%   Copyright (c) 2005 Gabriel Peyré


if nargin<3
    M = [];
end
if nargin<4
    plot_type = 1;
end

% be sure that there is a single QT
QT = QT(:,:,1);
Theta = Theta(:,:,1);

str = 'r';
str_geom = 'b';


hold on;

% display image
if ~isempty(M)
    imagesc([0 1],[0 1],M');
end
plot_square([0,0], 1, str);
axis square;
axis off;

j_min = min(QT(:));
j_max = max(QT(:));
n = size(QT,1);

% display subdivision
for j=log2(n):-1:j_min
    w = 2^j/n;
    for kx=0:n/2^j-1
        for ky=0:n/2^j-1
            pos = [kx,ky]*2^j/n;
            if QT(kx*2^j+1, ky*2^j+1)==j
               % this is a leaf
               if plot_type==1 || plot_type==2
                    plot_square_geometry( Theta(kx*2^j+1, ky*2^j+1), pos, w, str_geom );
               end
            elseif QT(kx*2^j+1, ky*2^j+1)<j
                % continue subdivision
                if plot_type==1 || plot_type==3
                    plot_cross(pos, w, str);
                end
            end
        end
    end
end


hold off;

axis tight;
axis square;
axis ij;
colormap gray(256); 


% CAUTION : since matlab image are X/Y swapper, 
% the display swap X and Y.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_cross(pos, w, str)

pos = swap_pos(pos);

if nargin<3
    str = 'r';
end

x = [pos(1), pos(1)+w];
y = [pos(2)+w/2, pos(2)+w/2];
plot(x,y, str);
x = [pos(1)+w/2, pos(1)+w/2];
y = [pos(2), pos(2)+w];
plot(x,y, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_square(pos, w, str)

% pos = swap_pos(pos);

if nargin<3
    str = 'r';
end

x = [pos(1), pos(1)+w, pos(1)+w, pos(1), pos(1)];
y = [pos(2), pos(2), pos(2)+w, pos(2)+w, pos(2)];
plot(x,y, str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_square_geometry(theta,pos,w, str)

if nargin<4
    str = 'b';
end

% pos = pos(2:-1:1);
x = pos(1)+w/2 + w/2*[cos(theta), -cos(theta)];
y = pos(2)+w/2 + w/2*[sin(theta), -sin(theta)];
plot(x,y, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos1 = swap_pos(pos)

pos1 = pos;
% pos1 = pos(2:-1:1);
%% pos1(1) = 1-pos1(1);