function [x_min, x_max, y_min, y_max] = mat2txt(filename, X, Y, Z)
%%
% this function is used to transfer the matlab surface matrix to the data
% that autoCAD (Solidworks) can read.  
%%

%% Saving direction

%filename = 'thickness_CAD';
saving_dir = '.\0.comsol_verification\';
file_dir = [saving_dir, filename];              % prefix of saving file name of retrieval phase distribution

%% Separate the matrix X, Y, Z into xx, yy, zz
% 
% % ============================== EXAMPLE ==================================
% % x = [1 2 3]=[x1 x2 x3];   y = [1 2 3]=[y1 y2 y3];
% % X = [1     2     3  = [x1  x2  x3      Y = [1     1     1 = [y1  y1  y1
% %      1     2     3     x1  x2  x3           2     2     2    y2  y2  y2
% %      1     2     3]    x1  x2  x3];         3     3     3]   y3  y3  y3];
% % Z = [(x1,y1) (x2,y1) (x3,y1)
% %      (x1,y2) (x2,y2) (x3,y2)
% %      (x1,y3) (x2,y3) (x3,y3)];
% % xx = [1  = [x1        yy = [1  = [y1         zz = [(x1,y1)
% %       1     x1              2     y2               (x1,y2)
% %       1     x1              3     y3               (x1,y3)
% %       2     x2              1     y1               (x2,y1)
% %       2     x2              2     y2               (x2,y2)
% %       2     x2              3     y3               (x2,y3)
% %       3     x3              1     y1               (x3,y1)
% %       3     x3              2     y2               (x3,y2)
% %       3]    x3];            3]    y3];             (x3,y3)];   
% % =========================================================================
% 
% [r, c] = size(z);
% x = linspace(1,c,c);
% y = linspace(1,r,r);
% [X, Y] = meshgrid(x(:), y(:));
% Z = z(:, :);
% 
% [r, c] = size(X);
% xx = reshape(X, r*c, 1);
% yy = reshape(Y, r*c, 1);
% zz = reshape(Z, r*c, 1);
% 
% %%
% 
% fid=fopen([file_dir,'.scr'],'w');
% fprintf(fid,'spline\n');
% for ii = 1:length(xx)
%     fprintf(fid,'%g,%g,%g\n',xx(ii),yy(ii),zz(ii));
% end
% fprintf(fid,'\n\n\nzoom\ne\n');
% fclose(fid);
% 

%% Separate the matrix X, Y, Z into vector xx, yy, zz

% X = [x1  x2  x3      Y = [y1  y1  y1      Z = [(x1,y1) (x2,y1) (x3,y1)
%      x1  x2  x3           y2  y2  y2           (x1,y2) (x2,y2) (x3,y2)
%      x1  x2  x3];         y3  y3  y3];         (x1,y3) (x2,y3) (x3,y3)];
%                                   ||
%                                   ||
%                                   ||
%                                   \/
%           xx = [x1        yy = [y1         zz = [(x1,y1)
%                 x2              y1               (x2,y1)
%                 x3              y1               (x3,y1)
%                 x1              y2               (x1,y2)
%                 x2              y2               (x2,y2)
%                 x3              y2               (x3,y2)
%                 x1              y3               (x1,y3)
%                 x2              y3               (x2,y3)
%                 x3];            y3];             (x3,y3)];   

[r, c] = size(X);
xx = reshape(X', r*c, 1);
yy = reshape(Y', r*c, 1);
zz = reshape(Z', r*c, 1);

%% Save in .txt file

fid = fopen([file_dir,'.txt'],'w');
fprintf(fid,'%%x y z\n');
for ii = 1:length(xx)
    fprintf(fid,'%g    %g    %g\n',xx(ii),yy(ii),zz(ii));
end
fclose(fid);

x_min = min(xx);
x_max = max(xx);
y_min = min(yy);
y_max = max(yy);


%%