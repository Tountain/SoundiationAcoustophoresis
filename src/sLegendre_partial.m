function [P_nm_partial_1, P_nm_partial_2, P_nm_partial_3] = sLegendre_partial(n,m,theta)
%%
% Obtain the deviation of associated Legendre function 
% n--lower order
% m--upper order
% theta--polar angle, radian mesurement (SHOULD be row vector)  [0,pi]
%%

P_nm_partial_1 = 0;
P_nm_partial_2 = 0;
P_nm_partial_3 = 0;
if m >= 0       
    
    P_nm = sLegendre(n,m,theta);                        % dependent variable is 'cos(theta)' not 'theta'!
    P_n1m = sLegendre(n+1,m,theta);                     % P_nm(cos(theta)), m is positive
    P_1nm = sLegendre(n-1,m,theta);
    P_nm1 = sLegendre(n,m+1,theta);
%     if m > n-1
%         P_1nm = 0;
%     else
%         P_1nm = sLegendre(n-1,m,theta);
%     end
%     if m+1 > n
%         P_nm1 = 0;
%     else
%         P_nm1 = sLegendre(n,m+1,theta);
%     end
    
    if length(theta) == 1  
        if (0.001/pi >= theta && theta >= 0)                % avoid the zero denominator '(1-cos(theta)^2)' in following calculations.
            theta = theta + 0.001/pi;
        elseif (theta >= pi-0.001/pi && theta <= pi)
            theta = theta - 0.001/pi;
        end
    else
        for ii = 1 : length(theta)
            if (0.001/pi >= theta(ii) && theta(ii) >= 0)    % avoid the zero denominator '(1-cos(theta)^2)' in following calculations.
                theta(ii) = theta(ii) + 0.001/pi;
            elseif (theta(ii) >= pi-0.001/pi && theta(ii) <= pi)
                theta(ii) = theta(ii) - 0.001/pi;
            end
        end
    end
    
    if length(theta) == 1                               % Eq.(6.36)
        P_nm_partial_1 = ((n+1)*cos(theta)*P_nm - (n-m+1)*P_n1m) / (1-cos(theta)^2);                       
    else
        P_nm_partial_1 = ((n+1)*cos(theta).*P_nm - (n-m+1)*P_n1m) ./ (1-cos(theta).^2);
    end
    if length(theta) == 1                               % Eq.(6.37)
        P_nm_partial_2 = (-n*cos(theta)*P_nm + (n+m)*P_1nm) / (1-cos(theta)^2);                       
    else
        P_nm_partial_2 = (-n*cos(theta).*P_nm - (n+m)*P_1nm) ./ (1-cos(theta).^2);
    end
    if length(theta) == 1                               % Eq.(6.38)
        P_nm_partial_3 = (-sqrt(1-cos(theta)^2)*P_nm1 - m*cos(theta)*P_nm) / (1-cos(theta)^2);                       
    else
        P_nm_partial_3 = (-sqrt(1-cos(theta).^2).*P_nm1 - m*cos(theta).*P_nm) ./ (1-cos(theta).^2);
    end

elseif m < 0    % Eq.(6.31) + Eq.(6.36) or Eq.(6.37) or Eq.(6.38) 
    
    m = abs(m);                              % transfer the negative m to positive m
    coeff = factorial(n-m)/factorial(n+m);
    P_nm = sLegendre(n,m,theta);             % P_n(-m)(cos(theta)), m is positive. dependent variable is 'cos(theta)' not 'theta'!
    P_n1m = sLegendre(n+1,m,theta);          % P_nm(cos(theta)), m is positive
    P_1nm = sLegendre(n-1,m,theta);
    P_nm1 = sLegendre(n,m+1,theta);
%     if m > n-1
%         P_1nm = 0;
%     else
%         P_1nm = sLegendre(n-1,m,theta);
%     end
%     if m+1 > n
%         P_nm1 = 0;
%     else
%         P_nm1 = sLegendre(n,m+1,theta);
%     end
    
    if length(theta) == 1  
        if (0.001/pi >= theta && theta >= 0)                % avoid the zero denominator '(1-cos(theta)^2)' in following calculations.
            theta = theta + 0.001/pi;
        elseif (theta >= pi-0.001/pi && theta <= pi)
            theta = theta - 0.001/pi;
        end
    else
        for ii = 1 : length(theta)
            if (0.001/pi >= theta(ii) && theta(ii) >= 0)    % avoid the zero denominator '(1-cos(theta)^2)' in following calculations.
                theta(ii) = theta(ii) + 0.001/pi;
            elseif (theta(ii) >= pi-0.001/pi && theta(ii) <= pi)
                theta(ii) = theta(ii) - 0.001/pi;
            end
        end
    end
    
    if length(theta) == 1                               % Eq.(6.31) + Eq.(6.36)
        P_nm_partial_1 = (-1)^(m) * coeff * ((n+1)*cos(theta)*P_nm - (n-m+1)*P_n1m) / (1-cos(theta)^2);                       
    else
        P_nm_partial_1 = (-1)^(m) * coeff * ((n+1)*cos(theta).*P_nm - (n-m+1)*P_n1m) ./ (1-cos(theta).^2);
    end
    if length(theta) == 1                               % Eq.(6.31) + Eq.(6.37)
        P_nm_partial_2 = (-1)^(m) * coeff * (-n*cos(theta)*P_nm + (n+m)*P_1nm) / (1-cos(theta)^2);                       
    else
        P_nm_partial_2 = (-1)^(m) * coeff * (-n*cos(theta).*P_nm - (n+m)*P_1nm) ./ (1-cos(theta).^2);
    end
    if length(theta) == 1                               % Eq.(6.31) + Eq.(6.38)
        P_nm_partial_3 = (-1)^(m) * coeff * (-sqrt(1-cos(theta)^2)*P_nm1 - m*cos(theta)*P_nm) / (1-cos(theta)^2);                       
    else
        P_nm_partial_3 = (-1)^(m) * coeff * (-sqrt(1-cos(theta).^2).*P_nm1 - m*cos(theta).*P_nm) ./ (1-cos(theta).^2);
    end
    
end


% reference equations can be found in Eq.(6.36), Eq.(6.37), Eq.(6.38)and
% Eq.(6.39) of E.G.W@1999@