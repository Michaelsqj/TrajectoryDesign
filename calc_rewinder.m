function [gr, kr] = calc_rewinder(k, g, gmax, smax, T)
    % This is the function to calculate the rewinder, using the maximum
    % gradient and slew rate
    % 
    % Input:
    %   kend: last kspace position, 3x1 (1/cm)
    %   gend: last gradient, 3x1
    %   gmax: G/cm
    %   smax: G/cm/ms
    %   Ts:   ms
    % Return:
    %   gr:   Nx3
    %   kr:   Nx3
    
    %---------------------------------------------
    %   Generate Rewiners 
    %     This is not the fasteest way to rewind gradients since it limits slew
    %     to gmax/sqrt(3) on each axis 
    %---------------------------------------------

    gmax3 = gmax/sqrt(3);
    smax3 = smax/sqrt(3);
    max_pts = 0;
%     readout_pts = size(k,1);
    max_rewindtime = 0;
    done = [0 0 0];
    for pass = 1:2
        for dir =1:3

            if done(dir)==1
                continue;
            end

            gX = g(:,dir);
            gx_ramp_width= abs(gX(end)) / smax3;
            gx_ramp_area = 0.5*gx_ramp_width*gX(end);
            areax_end = -sum(gX*T);

            disp(['Area end = ',num2str(areax_end),' gx_ramp_end=',num2str(gx_ramp_area)]);

            if  sign( areax_end - gx_ramp_area) == sign(gx_ramp_area)
                if pass==2
                    disp(['Dir=',num2str(dir),'Gradient end = ',num2str(gX(end)),' Area end = ',num2str(areax_end)]);
                    disp('Gradient should not be fully rewound - calculating bridged gradient');

                    w = 2*abs(areax_end/gX(end));
                    if w > max_rewindtime

                        % Calculate trap + square
                        w1 = 999999;
                        w2 = 999999;
                        w3 = 999999;
                        gend  = gX(end);
                        gmaxR =  gend + smax3*T*sign(gend);

                        while ( (w2+w1+w3) > max_rewindtime ) && ( abs(gmaxR) < gmax3);
                            disp(['W2 is ',num2str(w2+w),' of ',num2str(max_rewindtime)]);
                            gmaxR = gmaxR + smax3*T*sign(gend);

                            %Ramp down
                            w1  = abs(gmaxR / smax3);
                            pts = ceil( w1 / T);
                            gX_rampdn = linspace(gmaxR,0,pts+1);

                            %Ramp up
                            w2  = abs(gmaxR-gX(end)) / smax3;
                            pts = ceil( w2 / T);
                            gX_rampup = linspace(gX(end),gmaxR,pts+1);
                            %gX_rampup = gX_rampup(1:end);

                            %Flat top
                            area_target = T*( -sum(gX_rampup)-sum(gX)-sum(gX_rampdn));
                            w3 = area_target / gmaxR;
                            pts = ceil(w3/T);
                            gx_rewind_trap = ones(1,pts)*gmaxR;
                        end

                        g_rewind = [gX_rampup'; gx_rewind_trap'; gX_rampdn'];
                        area_adjust = areax_end - sum(T*g_rewind);
                        %Scale just the amount stronger than gX(end)
                        ramp_down = linspace(g_rewind(1),0,numel(g_rewind))';
                        gdiff = g_rewind - ramp_down;
                        gdiff = gdiff * ( areax_end - sum(ramp_down*T))/sum(T*gdiff);
                        g_rewind = gdiff + ramp_down;

                        %calculate constant an
                        gX = [gX; g_rewind];


                    else
                        disp('Ramp down only')
                        % Gradient is a ramp down only
                        pts = 1+ceil( w / T);
                        gX_rewind = linspace(1,0,pts);
                        a_rewind  = areax_end / (T*sum(gX_rewind));
                        gX = [gX; a_rewind*gX_rewind'];
                    end

                    g2{dir} = gX;
                    max_pts = max(numel(gX),max_pts);

                end

                disp(['Sum G = ',num2str(sum(T*gX))])
            elseif pass==1
                disp(['Dir=',num2str(dir),'Gradient end = ',num2str(g(end,dir)),' Area end = ',num2str(areax_end)]);
                disp(['Gradient Rewound + Trapezoid']);
                % Gradient maust be rewound first
                time_rewind = abs(g(end,dir)) / smax3;
                pts = ceil( time_rewind/ T);
                gX_rewind = linspace(g(end,dir),0,pts+1);
                gX_rewind = gX_rewind(2:end);
                gX = [gX; gX_rewind'];

                % Area
                area_trap =   -sum(gX*T);
                time_trap = sqrt( abs(area_trap) / smax3);
                pts = ceil( time_trap / T);
                gX_trap = [ linspace(0,1,pts) linspace(1,0,pts)];
                a_trap =  area_trap / (T*sum(gX_trap));
                if( abs( a_trap) > gmax3)
                    pts_up = ceil(gmax3/smax3/T);
                    pwa = pts_up*T;
                    pw = (abs(area_trap) - pwa*gmax3)/gmax3;
                    pts = ceil(pw/T);
                    gX_trap = [ linspace(0,1,pts_up) ones(1,pts) linspace(1,0,pts_up)];
                    a_trap =  area_trap / (T*sum(gX_trap));
                end

                gX = [gX; a_trap*gX_trap'];

                max_rewindtime = max( 2*time_trap+time_rewind, max_rewindtime);

                g2{dir} = gX;
                max_pts = max(numel(gX),max_pts);
                done(dir)=1;

                disp(['Sum G = ',num2str(sum(T*gX))])
            end

        end
    end

    % Append the Rewinder
    g3 = zeros(max_pts,3);
    g3(1:numel(g2{1}),1) = g2{1};
    g3(1:numel(g2{2}),2) = g2{2};
    g3(1:numel(g2{3}),3) = g2{3};
    g3 = [g3' zeros(3,1) ]'; % Add aditional dwell point
    g = g3;
    gr = g; 

    %Generate K-space from gradient
    kr = 4258*cumtrapz(g)*T/1000;

end