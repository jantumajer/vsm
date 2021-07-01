% Vaganov-Shashkin Tree Ring Growth Model (VSM) Parameters

% option to display program header and progress 
parameters.display = 1;         % set to 1 to write headers and program progress to screen, 0 for improved speed and no screen output

% Main Program and Growth Block Parameters 
% Piecewise growth function parameters for temperature
parameters.Tf(1) =      0.0000;  % minimum temperature (C) for growth                
parameters.Tf(2) =     15.0000;  % growth rate is max in the range T2-T3 (lower optimal temperature, C)   
parameters.Tf(3) =     22.0000;  % growth rate is max in the range T2-T3 (upper optimal temperature, C)
parameters.Tf(4) =     30.0000;  % maximum temperature (C) for growth

% Piecewise growth function parameters for soil moisture
parameters.Wf(1) =       .0200;   % minimum soil moisture for growth (v/v)
parameters.Wf(2) =       .2000;   % the growth rate is max in range W2-W3 (lower optimal soil moisture, v/v)
parameters.Wf(3) =       .8000;   % the growth rate is max in range W2-W3 (upper optimal soil moisture, v/v)
parameters.Wf(4) =       .9000;   % growth is stopped at this soil moisture (v/v)

parameters.SNo   =       .0000;  % initial snowpack (mm)
parameters.Wo    =       .2183;  % initial soil moisture (v/v)
parameters.rootd =    500.0000;   % the root/soil melt depth (mm)
parameters.rated =       .0200;  % the rate of water drainage from soil
parameters.Pmax  =     10.0000;  % maximum rate of infiltration water into soil  (mm/day)

% ///MODIFICATION///
parameters.Wmax  =       parameters.Wf(4);  % maximum soil moisture (field capacity, (v/v)
parameters.Wmin  =       parameters.Wf(1);  % minimum soil moisture (wilting point, v/v)

% interception and transpiration parameters
parameters.k(1)  =       .7000;  % k1 (1-k1) is the interception precipitation by the tree crown
parameters.k(2)  =       .1500;  % k2 coefficient for calculation the transpiration
parameters.k(3)  =       .0800;  % k3 coefficient for calculation the transpiration

% soil and snow melting parameters
parameters.Tm    =      20.000;  % sum of temperature for start soil melting (C)
parameters.a(1)  =     10.0000;  % a1, coefficient of soil melting
parameters.a(2)  =       .0060;  % a2, coefficient of soil melting
parameters.Tg    =     10.0000;  % sum of temperature to start growth (C)
parameters.SNr   =      5.0000;  % the rate of snow melting (mm/C/day)
parameters.SNmt  =       .0000;  % minimum temperature for snow melting

% Some switches and variable storage 
parameters.K(1)  =          0;  % soil melting switch: yes (1); no (0)
parameters.K(4)  =          1;  % snow melting switch: yes (1); no (0)
parameters.K(8)  =         50;  % Maximum duration (days) of latewood formation
parameters.K(9)  =         10;  % the period (days) over which to sum temperature to calculate start of growth
  % !!! Be aware that K(9) in this script equals tbeg-1 in VS-Oscilloscope

  parameters.K(10) =         10;  % the period (days) over which to sum temperature to calculate start soil melting

% Parameters for Cambial Block
parameters.ndc   =          1;  % Use previous cambium for following year? 1 = yes (dynamic cambium), 0 = no(static);

% Growth rate parameters
parameters.b(1)  =        .04;  % the critical growth rate (Vcr)
parameters.b(4)  =        0.6;  % %he correction of growth rate (Gr*b(4))
parameters.b(5)  =        1.9;  % The correction of Vo(j*b(5)) and Vmin(j*b(5))
parameters.b(6)  =        .25;  % Vo(j)= b(6)*j+b(7)
parameters.b(7)  =        .42;  % b(7) and b(6) determinate the fuction  Vo(j)
parameters.b(8)  =       2.50;  % Vmin(j)=(EXP((b(8)+j)*0.4)-5.0)*b(9)
parameters.b(9)  =       0.04;  % b(8) and b(9) determine the fuction Vmin(j)
parameters.b(10) =       1.00;  % The growth rate in the S, G2 & M phases 
parameters.b(11) =       8.00;  % The maximum size of a cell in the G1 phase (SIG1)
parameters.b(12) =       9.00;  % The maximum size of a cell in the S phase
parameters.b(13) =       9.50;  % The maximum size of a cell in the G2 phase
parameters.b(14) =      10.00;  % The maximum size of a cell in the M phase (SIM)
parameters.b(15) =        1.0;  % The time step in cambium model (1/b(15)) = number of c-cycles/day)

% ///MODIFICATION/// - new input parameters 
parameters.Vs      =      .04; % The critical growth rate for environmental block - cambial activity does not ceasse if Gr>Vs in two consecutive days in a row even if cumulative temperature drops below Tg/K(9) threshold
parameters.cambial =        1; % Binary variable - should cambial block be run (1) or not (0=TRW chronology simulated from annual sums of GR as in VS-Oscilloscope)