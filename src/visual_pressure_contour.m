function [figure_num] = visual_pressure_contour(gr, pt_v, pi_v, ps_v, prs_v, xx, yy, zz)
%%
% this function is used to visualize the pressure data under different
% position of Cartesian coordinate.
%%

parameters;

% pi_v = pi_v / max(abs(pi_v));

%% visualize in xz-plain

if abs(min(yy)-max(yy)) < (particle_radius/4)      % in xz-plain, y == 0

    figure_num = 3;
    figure(figure_num);

    X1 = linspace(min(xx),max(xx),gr);
    Z1 = linspace(min(zz),max(zz),gr);
    [XXI, ZZI, PI_V] = griddata(xx(1:10:end),zz(1:10:end),pi_v(1:10:end),X1',Z1,'v4');
    [XXS, ZZS, PS_V] = griddata(xx(1:10:end),zz(1:10:end),ps_v(1:10:end),X1',Z1,'v4');
    [XXRS, ZZRS, PRS_V] = griddata(xx(1:10:end),zz(1:10:end),prs_v(1:10:end),X1',Z1,'v4');
%     X1 = linspace(min(xx),max(xx),gr*5);
%     Z1 = linspace(min(zz),max(zz),gr*5);    
%     [XXI, ZZI, PI_V] = griddata(xx(1:1:end),zz(1:1:end),pi_v(1:1:end),X1',Z1,'cubic');
%     [XXS, ZZS, PS_V] = griddata(xx(1:1:end),zz(1:1:end),ps_v(1:1:end),X1',Z1,'cubic');
%     [XXRS, ZZRS, PRS_V] = griddata(xx(1:1:end),zz(1:1:end),prs_v(1:1:end),X1',Z1,'cubic');
    XXT = (XXI + XXS + XXRS)/3;            % avoid the numerical "eating"
    ZZT = (ZZI + ZZS + ZZRS)/3;
    PT_V = PI_V + PS_V + PRS_V;             
    
    %contourf(XXT, ZZT, PT_V);
    pclr = pcolor(real(XXI), real(ZZI), real(PI_V));
    set(pclr, 'LineStyle','none');
    %set(pclr, 'box', 'on');
    colorbar;
    caxis([-300,300]);

    xlabel('X [m]');
    ylabel('Z [m]');
    
    global N
    title(['Truncation Number N = ', num2str(N)]);
    
end


%% visualize in xy-plain

if abs(min(zz)-max(zz)) < (particle_radius/4)      % in xy-plain, z == fixed_z*pa

    figure_num = 2;
    figure(figure_num);

    X1 = linspace(min(xx),max(xx),gr);
    Y1 = linspace(min(yy),max(yy),gr);    
    [XXI, YYI, PI_V] = griddata(xx(1:10:end),yy(1:10:end),pi_v(1:10:end),X1',Y1,'v4');
    [XXS, YYS, PS_V] = griddata(xx(1:10:end),yy(1:10:end),ps_v(1:10:end),X1',Y1,'v4');
    [XXRS, YYRS, PRS_V] = griddata(xx(1:10:end),yy(1:10:end),prs_v(1:10:end),X1',Y1,'v4');
%     X1 = linspace(min(xx),max(xx),gr*5);
%     Y1 = linspace(min(yy),max(yy),gr*5);
%     [XXI, YYI, PI_V] = griddata(xx(1:1:end),zz(1:1:end),pi_v(1:1:end),X1',Y1,'cubic');
%     [XXS, YYS, PS_V] = griddata(xx(1:1:end),zz(1:1:end),ps_v(1:1:end),X1',Y1,'cubic');
%     [XXRS, YYRS, PRS_V] = griddata(xx(1:1:end),zz(1:1:end),prs_v(1:1:end),X1',Y1,'cubic');
    XXT = (XXI + XXS + XXRS)/3;            % avoid the numerical "eating"
    YYT = (YYI + YYS + YYRS)/3;
    PT_V = PI_V + PS_V + PRS_V;             

    %contourf(XXT, ZZT, PT_V);
    pclr = pcolor(real(XXT), real(YYT), real(PT_V));
    set(pclr, 'LineStyle','none');
    %set(pclr, 'box', 'on');
    colorbar;
    %caxis([-0.6,1.2]);

    xlabel('X [m]');
    ylabel('Y [m]');
    title(['XY-plain on Z = ', num2str(round(mean(zz)*10^6)), ' \mum']);

end
 
%% saving XXT, YYT (or ZZT), PI_V, PS_V, PRS_V and PT_V for Plane wave
%% expansion (PWE) method, inactive if the PWE is not used

% if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer')) == 1
%     ii = 1;
%     while 1
%         dir = ['.\data_pres_transducer_', fluid, '\', num2str(ii)];
%         if exist([dir, '.mat']) ~= 0        % if already exist the database, 
%              ii = ii + 1;                   % then do not create again for saving time.
%         else
%             break;
%         end
%     end
%     if abs(min(yy)-max(yy)) < (particle_radius/4)      % in xz-plain, y == 0
%         save([dir, '.mat'], 'XXT', 'ZZT', 'PI_V', 'PS_V', 'PRS_V', 'PT_V');
%     end
%     if abs(min(zz)-max(zz)) < (particle_radius/4)      % in xy-plain, z == fixed_z*pa
%         save([dir, '.mat'], 'XXT', 'YYT', 'PI_V', 'PS_V', 'PRS_V', 'PT_V');
%     end
% end

%%
