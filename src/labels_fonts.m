figure_FontSize = 16;

%% 3D

% title(['\it{\fontname{Times new roman}E_{\rm{RMS}}} \rm{\fontname{Times new roman} = } ', num2str(roundn(Err,-3))], ...
%     'FontSize', figure_FontSize);
% xlabel('\rm{\fontname{Times new roman}\it{x}-axis}', 'rotation', 22);
% ylabel('\rm{\fontname{Times new roman}\it{y}-axis}', 'rotation', -20);
% zlabel('\rm{\fontname{Times new roman}\it{ G(x,y)}}');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-8);

%% 2D

% title(['\it{\fontname{Times new roman}E_{\rm{RMS}}} \rm{\fontname{Times new roman} = } ', num2str(roundn(Err,-3))], ...
%     'FontSize', figure_FontSize);
xlabel('\rm{\fontname{Times new roman}\it{x}\rm{-axis [mm]}}');
ylabel('\rm{\fontname{Times new roman}\it{y}\rm{-axis [mm]}}');
set(get(gca,'XLabel'),'FontSize',figure_FontSize+8);
set(get(gca,'YLabel'),'FontSize',figure_FontSize+8);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize);

%% ylabel adjusting

set(gca,'XTick', [-0.02, -0.01, 0, 0.01, 0.02]);
set(gca,'XTickLabel', [-20, -10, 0, 10, 20]);
set(gca,'YTick', [-0.02, -0.01, 0, 0.01, 0.02]);
set(gca,'YTickLabel', [-20, -10, 0, 10, 20]);

%% colorbar

h1 = colorbar;                       % 一定要放置于 hc = findobj(allchild(gcf), 'Type', 'axes'); hc2 = findobj(allchild(hc), 'Type', 'text');之后。
% colormap(color);                   % 因为colorbar得加入，将使得'axes'得'text'对象得第三个和第四个不再仅仅代表X轴和Y
% set(h1,'YTick',[1 3 5 7 9 11]);
% set(h1,'YTickLabel',{'0','10','20','30','40', '50'});
set(get(h1,'title'),'string','\rm{\fontname{Times new roman}{Pa}}','fontsize',figure_FontSize+1);
set(findobj('FontSize',10),'FontSize',figure_FontSize);


%%