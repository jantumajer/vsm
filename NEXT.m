%%% Calculation of expected growth rate for the next day (t+1)
% Estimates the growth rate of the next day based on known temperature and daylength conditions and unknown moisture level
% Assumes stable moisture: sm(t, cyear) == sm(t+1), cyear)

function[Grnext] = NEXT(t,T,P,sm,ndl,cyear,ndays,parameters)
        % First, I need GrE time series
            if ndays(cyear)==366
                GrE=ndl(:,2);
            else
                GrE=ndl(:,1);
            end

        % Then, I calculate growth rates to temperature in day t and t+1
        if t+1 > ndays(cyear)
            Grnext = 0;
          else
            
        x = T(t+1,cyear); y = T(t,cyear);
        if (x < parameters.Tf(1))
            GrT1(t,cyear) = 0;
        elseif (x >= parameters.Tf(1)) && (x <= parameters.Tf(2))
            GrT1(t,cyear) = (x - parameters.Tf(1))/(parameters.Tf(2) - parameters.Tf(1));
        elseif (x >= parameters.Tf(2)) && (x <= parameters.Tf(3))
            GrT1(t,cyear) = 1;
        elseif (x >= parameters.Tf(3)) && (x <= parameters.Tf(4))
            GrT1(t,cyear) = (parameters.Tf(4) - x)/(parameters.Tf(4) - parameters.Tf(3));
        else
            GrT1(t,cyear) = 0;
        end
        
        if (y < parameters.Tf(1))
            GrT(t,cyear) = 0;
        elseif (y >= parameters.Tf(1)) && (y <= parameters.Tf(2))
            GrT(t,cyear) = (y - parameters.Tf(1))/(parameters.Tf(2) - parameters.Tf(1));
        elseif (y >= parameters.Tf(2)) && (y <= parameters.Tf(3))
            GrT(t,cyear) = 1;
        elseif (y >= parameters.Tf(3)) && (y <= parameters.Tf(4))
            GrT(t,cyear) = (parameters.Tf(4) - y)/(parameters.Tf(4) - parameters.Tf(3));
        else
            GrT(t,cyear) = 0;
        end
        
        % Next, growth rate to soil moisture in day t
        x = sm(t,cyear);
        if (x < parameters.Wf(1))
            GrW(t,cyear) = 0;
        elseif (x >= parameters.Wf(1)) && (x <= parameters.Wf(2))
            GrW(t,cyear) = (x - parameters.Wf(1))/(parameters.Wf(2) - parameters.Wf(1));
        elseif (x >= parameters.Wf(2)) && (x <= parameters.Wf(3))
            GrW(t,cyear) = 1;
        elseif (x >= parameters.Wf(3)) && (x <= parameters.Wf(4))
            GrW(t,cyear) = (parameters.Wf(4) - x)/(parameters.Wf(4) - parameters.Wf(3));
        else
            GrW(t,cyear) = 0;
        end
        
        %%%% Soil moisture model component to estimate soil moisture level for the next day %%%%
        % transpiration
        trt             = parameters.k(2) * exp(T(t,cyear) * parameters.k(3));
        trans(t,cyear)  = trt * (GrE(t) * min(GrT(t,cyear),GrW(t,cyear)));
        
        prrain(t, cyear) = P(t,cyear);

        % calculate soil moisture for the next day
        % because of late summer-autumn season of cambial activity cessation I drop snow-melt
        K(6) = min(parameters.k(1) * prrain(t, cyear), parameters.Pmax);
        w    = sm(t,cyear) * parameters.rootd * (1 - parameters.rated) + K(6) - trans(t,cyear);  % available water used for soil moisture calculation
        % w    = sm(t,cyear) * rootd * (1 - rated) + K(6) - trans + xm;  % available water used for soil moisture calculation
        if isnan(w)==1; w = 0; end            % error catching
        if w <= 0; w = 0; end                 % available snow can't be less than zero
        
        sm(t+1,cyear) = w/parameters.rootd;
                if sm(t+1,cyear) <= parameters.Wmin; sm(t+1,cyear) = parameters.Wmin; end;   % error catching
                if sm(t+1,cyear) >= parameters.Wmax; sm(t+1,cyear) = parameters.Wmax; end;   % error catching
                if isnan(sm(t+1,cyear))==1; sm(t+1,cyear) = parameters.Wmin; end; % error catching
        %%%%   end of soil moisture model %%%%
        
        % Then I can calculate growth rate to soil moisture in day t+1
        y = sm(t+1,cyear);
        if (y < parameters.Wf(1))
            GrW1(t,cyear) = 0;
        elseif (y >= parameters.Wf(1)) && (y <= parameters.Wf(2))
            GrW1(t,cyear) = (y - parameters.Wf(1))/(parameters.Wf(2) - parameters.Wf(1));
        elseif (y >= parameters.Wf(2)) && (y <= parameters.Wf(3))
            GrW1(t,cyear) = 1;
        elseif (y >= parameters.Wf(3)) && (y <= parameters.Wf(4))
            GrW1(t,cyear) = (parameters.Wf(4) - y)/(parameters.Wf(4) - parameters.Wf(3));
        else
            GrW1(t,cyear) = 0;
        end
        
        
        % daily growth rate calculation
        Grnext = GrE(t+1) * min(GrT1(t,cyear),GrW1(t,cyear));
  endif