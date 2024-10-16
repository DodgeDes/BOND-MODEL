% This is the final model for the implementation of the bond number. 
% The model was simulated for at 2D reservoirs with and without the effect
% of the bond number. Using a constant permeability and critical condensate
% saturation of 0.2:
% This model is to simulated different critical condensate saturation
% (0.05, 0.1, 0.3, 0.4) and also the base model was 0.2. The value of 0.5
% was used for Tau in all the simulations
%

clc; clear;
mrstModule add compositional ad-core ad-props mrst-gui linearsolvers; 
tic
% Input positive numbers
% nums = [1,2,3,4,5,10,13,16,18,21,31,70,130,300];
nums = [0.5,1.0,1.5,2,2.5,3,3.5,4,4.5,5,5.5,10,15,20,25,30,35,40,50,80,100,500];
% Initialize an empty array to store results
result = [];

% Loop through each number in the input array
for i = 1:length(nums)
    % Add the negative counterpart of the current number
    result = [result, -nums(i)];
end

% Append the original positive numbers to the result array
result = sort([result, nums]);

% Display the result
% disp(result);

x = result;
y = result;

% Create the tensor grid and plot it
G = tensorGrid(x, y);
G = computeGeometry(G);
% plotGrid(G); 
% % axis([-70 70 -70 70]);
% axis([-500 500 -500 500]);

%% DEFINE THE PETROPHYSICAL PROPERTIES
perm = 100; %in mD
% Define initial relative permeability values
kr_initial = perm * milli*darcy * ones(G.cells.num, 1); % Initial relative permeability for all grids

% Update the rock properties with initial relative permeability values
rock = makeRock(G, kr_initial, 0.2);
pv = poreVolume(G, rock);

%% SET UP THE COMPOSITIANL FLUID 
% Set up compositional fluid model
% Three-component system of Methane, CO2 and n-Decane
% names = {'Methane'  'n-Butane'  'n-Decane'};
names = {'CO2'  'C1'  'C2-C3' 'C4-C6'  'C9'  'C22'};

% Critical temperatures # Kelvin
% Tcrit = [190.6 425.2 622.1];
Tcrit = [304.2 190.6 330.5 453.8 606 869.7];
% Critical pressure 
% Pcrit = [4600155 3799688 2534138.299];
Pcrit = [7376460 4600155 4630552.5 3536242.5 2634450 1509742.5 ];
% Critical volumes
% Vcrit = [9.8628e-05 2.55e-04 6.0976e-04];
Vcrit = [0.000094 0.000099 0.00016964 0.0002986 0.0004904 0.0010589];
% Acentric factors
% acentricFactors = [0.008 0.193 0.443774];
acentricFactors = [ 0.225 0.008 0.119 0.226 0.359 0.788];
% Mass of components (as MRST uses strict SI, this is in kg/mol)
% molarMass = [0.016043 0.058124 0.134];
molarMass = [0.04401 0.016043 0.036 0.072 0.128 0.31];
% Initialize fluid
fluid = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass);

% bic = [ 0	0.01474853	0.0443719; 0.01474853	0	0.008452228; 0.0443719	0.008452228	0];
bic = [ 0	0 0  0  0 0; 0	0 0  0  0 0;0	0 0  0  0 0;0	0 0  0  0 0;0	0 0  0  0 0;0	0 0  0  0 0];
zcomp = [0.1618000 0.7098000 0.0790000 0.0260000 0.0200000 0.0034000];

fluid = setBinaryInteraction(fluid,bic);

%% SET UP INITIAL FLUID CONDITIONS
% This cell group flashes initial fluid composition to surface conditions
eosname = 'prcorr';
EOSModel = EquationOfStateModel(G, fluid, eosname);

%Surface Conditions for computation of density
p_sc = 100000*Pascal;  %atmospheric pressure
T_sc = (273.15 + 20)*Kelvin; % 20 celcius
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, zcomp, EOSModel);

%% SET UP THE BOUNDARY CONDITIONS
% outletP = 19*10^6*Pascal;
outletP = 22*10^6*Pascal;
inletP = 29*10^6*Pascal;
reservoirP = inletP;

bc = [];
bc = pside(bc, G, 'xmin', inletP, 'sat', [0, 1]);
bc = pside(bc, G, 'ymin', inletP, 'sat', [0, 1]);
bc = pside(bc, G, 'xmax', inletP, 'sat', [0, 1]);
bc = pside(bc, G, 'ymax', inletP, 'sat', [0, 1]);
bc.face = unique(bc.face); 
bc.components = repmat([0.1618000 0.7098000 0.0790000 0.0260000 0.0200000 0.0034000], numel(bc.face), 1);
% 
% bc = fluxside(bc, G, 'xmin', 1, 'sat', [0, 1]);  % 0 is oil and 1 is gas. confirmed from how the results are outputted
% bc = fluxside(bc, G, 'xmax', 1, 'sat', [0, 1]);
% bc = fluxside(bc, G, 'ymin', 1, 'sat', [0, 1]);
% bc = fluxside(bc, G, 'ymax', 1, 'sat', [0, 1]);
% bc.face = unique(bc.face); %this was commented now to see if there will be
% % bc.components = repmat([0.8325,0.1125, 0.055], numel(bc.face), 1);
% bc.components = repmat([0.1618000 0.7098000 0.0790000 0.0260000 0.0200000 0.0034000], numel(bc.face), 1);

% cellsWell = ceil(prod(dims) / 2);
cellsWell = round(G.cells.num / 2);
W = [];
well_rate = -0.35; %m3/s


% W = addWell(W,G, rock, 100, 'Type', 'bhp','comp_i', [0.5,0.5], 'Name', 'Pro', 'Val', outletP,'radius',0.0762, 'sign', -1);this was commented now to see if there will be
% W = addWell(W,G, rock, 100, 'Type', 'bhp','comp_i', [0.5,0.5], 'Name', 'Pro', 'Val', outletP, 'sign', -1);
W = addWell(W,G, rock, cellsWell, 'Type', 'bhp','comp_i', [0.5,0.5], 'Name', 'Pro', 'Val', outletP, 'sign', -1);
% W = addWell(W,G, rock, cellsWell, 'Type', 'rate','comp_i', [0.5,0.5], 'Name', 'Pro', 'Val', well_rate, 'sign', -1);


W(1).components = zcomp;


%% SET UP THE INITIAL FLUID MODEL
% this was commented to see if there will be a change.
% flowfluid = initSimpleADIFluid( ...
%                         'phases', 'OG', ...
%                         'mu',[1,1]*centi*poise, ...  % viscosity is just an input. It will be overwritten by the initialisation model
%                         'rho',[rhoO_S,rhoG_S], ...  % same as viscosity
%                         'n', [2,2], ...
%                         'c',[0,0]);     

%% SET UP THE INITIAL FLUID MODEL
flowfluid = initSimpleADIFluid( ...
                        'phases', 'OG', ...
                        'mu',[0.0758,0.0758]*centi*poise, ...  % viscosity is just an input. It will be overwritten by the initialisation model
                        'rho',[rhoO_S,rhoG_S]);      % same as viscosity

%% MAKE EACH CELL AND INDIVIDUAL REGION AND ASSIGN RELATIVE PERM TO EACH CELL
numCells = G.cells.num;
% Each cell is its own region
reg = (1:G.cells.num)';

% Update the fluid properties based on the random numbers
flowfluid.krO = cell(numCells, 1);
flowfluid.krG = cell(numCells, 1);

% OLD
% constants for relative permeability functions.
no = 2;
ng = 2;
kor = 0.2;     % residual oil sat
kgr = 0;       % residual gas sat
kroMax = 1.0;  % Maximum sat of oil
krgMax = 1.0;  % Maximum sat of gas

no_init = no;
ng_init = ng;
kor_init = kor;        % Residual oil saturation
kgr_init = kgr;        % Residual gas saturation
kroMax_init = kroMax;  % Maximum saturation of oil
krgMax_init = krgMax;  % Maximum saturation of gas

krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor+ kgr);
krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);

for i = 1:numCells
    % Assign the relative permeability function
    flowfluid.krO{i} = krOG;
    flowfluid.krG{i} = krG;
end


%% SET UP THE MODEL
arg = {G,rock,flowfluid,fluid, 'water', false};
sparse_backend = SparseAutoDiffBackend();
constructor = @GenericOverallCompositionModel;
modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelSparseAD = modelSparseAD.validateModel();
kr = modelSparseAD.FlowPropertyFunctions.RelativePermeability;
kr = kr.setRegions(reg);
modelSparseAD.FlowPropertyFunctions.RelativePermeability = kr;

%% SET THE INITIAL CONDDITION OF THE MODEL
ncomp = fluid.getNumberOfComponents();
s0 = [1, 1];
% T = 273.15 + 61.85*Kelvin;
T = 273.15 + 30*Kelvin;
state0 = initCompositionalState(G, reservoirP, T,s0,zcomp, modelSparseAD.EOSModel);
state0.T = repmat(T, G.cells.num, 1);
state0.components = repmat(zcomp, G.cells.num, 1);

%% First alternative: single separator

s = EOSSeparator('pressure', 100000*Pascal, 'T', 273.15+22*Kelvin);
sg = SeparatorGroup(s);
sg.mode = 'moles';
modelSparseAD.FacilityModel.SeparatorGroup = sg;

linsolve = selectLinearSolverAD(modelSparseAD);
disp(linsolve)


%% SET UP SCHEDULE
% num_iterations =50;  % In days
% tstep = 1; % for plot
% step =  tstep*day;  % for iteration
% totTime = num_iterations * day;  % total simulation time for 12 days
% dt = rampupTimesteps(totTime, step, 1);
num_iterations = 200;  % In days
tstep = 1; % for plot
step =  tstep/4*day;  % for iteration
totTime = num_iterations * day;  
dt = rampupTimesteps(totTime, step, 1);
% schedule = simpleSchedule(dt, 'W', W);
schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
% schedule = simpleSchedule(dt, 'bc', bc);
nls            = NonLinearSolver('LinearSolver', linsolve);
linsolve_block = AMGCL_CPRSolverBlockAD('tolerance', 1e-4);

gravity reset on

%% INITIALIZE DATA TO STORE COMPUTED DATA
own_states = cell(num_iterations, 1);
own_states{1} = state0;
own_wellSols = cell(num_iterations, 1);


% %% INPUT DATA FOR COMPUTATION FOR 3 COMPONENT FLUID
% pressureData = [35, 33,31,29,27,25,23,21,19,17,15]*1e6;  % PASCAL
% IFTData = [0,0,0,0.0006,0.0204,0.0792,0.1973,0.4047,0.7443,1.2753,2.0716]*0.001; % N/m (the issue is the ift data, check it and update it)
% 
%% IFT FOR SPE FLUID
% pressureData = [35, 33,31,29,27,25,23,21,19,17,15]*1e6;  % PASCAL
% IFTData = [0,0,0,0.0006,0.0204,0.0792,0.1973,0.4047,0.7443,1.2753,2.0716]*0.001; % N/m (the issue is the ift data, check it and update it)


% pressureData = [32, 30,29,28,26,24,22,20,18,16,14]*1e6;  % PASCAL
% IFTData = [0.00E+00,1.34E-04,1.65E-04,2.01E-04,2.90E-04,4.07E-04,5.63E-04,7.71E-04,1.06E-03,1.49E-03,2.12E-03]*0.001; % N/m (the issue is the ift data, check it and update it)

pressureData = [32, 31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14]*1e6;
IFTData = [0.00E+00,0.00E+00,0.00E+00,0.00E+00,4.53E-05,6.03E-05,7.86E-05,1.01E-04,1.27E-04,1.58E-04,1.96E-04,2.41E-04,2.95E-04,3.62E-04,4.49E-04,5.66E-04,7.29E-04,9.52E-04,1.25E-03];


% CN_base = 8.62*10^(-12);
BN_base = 3.0*10^(-6);
% BN_base = 10000000;
CN_base = 10000000;
tau = 4;

Bond_all = cell(num_iterations, 1);
Capillary_all = cell(num_iterations, 1);
IFT_C = cell(num_iterations, 1);
N_BN_all = cell(num_iterations, 1);
N_CN_all = cell(num_iterations, 1);
krr_all = cell(num_iterations, 1);
drho_all = cell(num_iterations, 1);
rp_type_all = cell(num_iterations, 1);
VMag_all = cell(num_iterations, 1);


%% SOLVER
for i = 1:num_iterations
    % Determine the timesteps corresponding to the current day
    start_idx = (i - 1) * tstep + 1;
    end_idx = i * tstep;

    % Modify the schedule to only contain the current day's timesteps
    single_step_schedule = schedule;
    single_step_schedule.step.val = schedule.step.val(start_idx:end_idx);
    single_step_schedule.step.control = schedule.step.control(start_idx:end_idx);

    % Simulate the schedule
   [wellSols, states, report] = simulateScheduleAD(state0, modelSparseAD, single_step_schedule, 'nonlinearsolver', nls);


    %% Solve for bond number and update the relative permeability
     % Interpolate IFT function
    cell_pressure = states{end}.pressure;  
    IFTInterpFunc = @(pressure) interp1(pressureData, IFTData, pressure, 'makima', 'extrap');
    ift_c = IFTInterpFunc(cell_pressure);


    %% Compute the Bond number
    mDtom = 9.869233*10^(-16);
    % L = 1.0e-9;
    L= perm*mDtom;
    g = 9.8;
    epsilon = 1.0e-30;

    rho_g = states{end}.PVTProps.Density{2};    % GAS
    rho_o = states{end}.PVTProps.Density{1};    % OIL
    drho = abs(rho_o - rho_g);

    Bond = zeros(numCells, 1);
    N_BN = zeros(numCells, 1);
    N_CN = zeros(numCells, 1);
    krr = zeros(numCells, 1);
    rp_type = zeros(numCells, 1);
    Capillary = zeros(numCells, 1);

    %% Compute the capillary number
    % Compute velocity
    v = faceFlux2cellVelocity(G, states{end}.flux(:, 1));
    v_mag = sqrt(sum(v.^2, 2));                 %velocity magniture
    mu = states{end}.PVTProps.Viscosity{1};     % Compute the viscosity of the gas phase{1}

    % % Compute capillary number
    % Ca = ( mu .* v_mag) ./ ift_c;   %internally computed capillary number
    % Nc = Nc_base ./ Ca;             % the normalized capillary number. If less than 1 then use the effect of capillary number. If more than 1 then used the base case RP.

    % Preallocate Bond array for each cell
    numCells = G.cells.num;
    % Each cell is its own region
    reg = (1:G.cells.num)';


    % Update the fluid properties based on the random numbers
    flowfluid.krO = cell(numCells, 1);
    flowfluid.krG = cell(numCells, 1);

    
    % if i > 1
    for j = 1:numCells
        % Compute Bond Number and capillary number
        if drho(j) > 0    % This takes care of the regions with single phase flows

            Bond(j) = (drho(j) * g * L) / (ift_c(j)+epsilon);
            Capillary(j) = (mu(j) * v(j)^2)/ (ift_c(j)+epsilon);

            N_BN(j) = BN_base / Bond(j);
            N_CN(j) = CN_base / Capillary(j);

            current_Nc =  N_CN(j);
            current_Nb =  N_BN(j);

            % Check the Bond number condition
            % if Bond(j) <= 0 || isnan(Bond(j))
            % if Bond(j) == 0 || isnan(Bond(j))
            if N_BN(j) < 1  % for the region where bond number will take effect

                mp = 1/tau;
                Xp(j) = 1 - exp(-mp * N_BN(j));
                % Xp(j) = 1 - exp(-mp);
                Xp_current = Xp(j);

                kor = kor_init * Xp(j);
  

                % kor = kor_init - kor_init * exp(-tau * N_BN(j));

                % Constants for relative permeability functions
                % no = no_init;
                % ng = ng_init;
                % Constants for relative permeability functions
                no = no_init;
                ng = ng_init;
                kgr = kgr_init;    % Residual gas saturation
                kroMax = kroMax_init;  % Maximum saturation of oil
                krgMax = krgMax_init;  % Maximum saturation of gas
                sr_tot = kor + kgr;

                krOG = coreyPhaseRelpermAD_DesNb(no, kor, kroMax, sr_tot,current_Nb,Xp_current) ;
                krG =  coreyPhaseRelpermAD_DesNb(ng, kgr, krgMax, sr_tot,current_Nb,Xp_current) ;

                % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
                % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);

                % % Constants for relative permeability functions
                % no = no_init;
                % ng = ng_init;
                % kor = kor_init;        % Residual oil saturation
                % kgr = kgr_init;        % Residual gas saturation
                % kroMax = kroMax_init;  % Maximum saturation of oil
                % krgMax = krgMax_init;  % Maximum saturation of gas
                % 
                % krOG = coreyPhaseRelpermAD_DesNc(no, kor, kroMax, kor+ kgr, current_Nb);
                % krG = coreyPhaseRelpermAD_DesNc(ng, kgr , krgMax, kgr+kor, current_Nb);

                rp_type(j) = 1;

            elseif N_CN(j) < 1  % for the region where capillary number will take effect 

                % Constants for relative permeability functions
                no = no_init;
                ng = ng_init;
                kor = kor_init;        % Residual oil saturation
                kgr = kgr_init;        % Residual gas saturation
                kroMax = kroMax_init;  % Maximum saturation of oil
                krgMax = krgMax_init;  % Maximum saturation of gas

                krOG = coreyPhaseRelpermAD_DesNc(no, kor, kroMax, kor+ kgr, current_Nc);
                krG = coreyPhaseRelpermAD_DesNc(ng, kgr , krgMax, kgr+kor, current_Nc);

                rp_type(j) = 2;

            else
                no = no_init;
                ng = ng_init;
                kor = kor_init;        % Residual oil saturation
                kgr = kgr_init;        % Residual gas saturation
                kroMax = kroMax_init;  % Maximum saturation of oil
                krgMax = krgMax_init;  % Maximum saturation of gas


                % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
                % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr);
                % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
                % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);
                krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
                krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);

                 rp_type(j) = 0;
            end

        else
            % for single phase flow
            no = no_init;
            ng = ng_init;
            kor = kor_init;        % Residual oil saturation
            kgr = kgr_init;        % Residual gas saturation
            kroMax = kroMax_init;  % Maximum saturation of oil
            krgMax = krgMax_init;  % Maximum saturation of gas


            % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
            % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr);
            % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
            % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);
            krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
            krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);
            rp_type(j) = 0;
        end 
        krr(j) = kor;

        % rp_type_collector = [];
        % rp_type_collector = [rp_type_collector, rp_type]';


        %compute the realtive permebility functions
        % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor + kgr);
        % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);

        % Assign the relative permeability function to the current cell
        flowfluid.krO{j} = krOG;
        flowfluid.krG{j} = krG;
    end
    % else
        % % for single phase flow
        % no = no_init;
        % ng = ng_init;
        % kor = kor_init;        % Residual oil saturation
        % kgr = kgr_init;        % Residual gas saturation
        % kroMax = kroMax_init;  % Maximum saturation of oil
        % krgMax = krgMax_init;  % Maximum saturation of gas
        % 
        % 
        % krOG = coreyPhaseRelpermAD_Des(no, kor, kroMax, kor+ kgr);
        % krG = coreyPhaseRelpermAD_Des(ng, kgr, krgMax, kgr+kor);
        % 
        % for k = 1:numCells
        %     % Assign the relative permeability function
        %     flowfluid.krO{k} = krOG;
        %     flowfluid.krG{k} = krG;
        % end
    % end

    %% store the well and grid information computed from this time step
    own_wellSols{i} = states{end}.wellSol;
    own_states{i} = states{end};

    arg = {G,rock,flowfluid,fluid, 'water', false};
    sparse_backend = SparseAutoDiffBackend();
    constructor = @GenericOverallCompositionModel;
    modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
    modelSparseAD = modelSparseAD.validateModel();
    kr = modelSparseAD.FlowPropertyFunctions.RelativePermeability;
    kr = kr.setRegions(reg);
    modelSparseAD.FlowPropertyFunctions.RelativePermeability = kr;


    %% Output computed values of the bond number and relative permeabilty
    Bond_all{i} = Bond;  % Store the Bond numbers for all cells in the current iteration
    Capillary_all{i} = Capillary;
    IFT_C{i} = ift_c;  % Store the IFT values for all cells in the current iteration
    N_BN_all{i} = N_BN;
    N_CN_all{i} = N_CN;
    krr_all{i} = krr;
    drho_all{i} = drho;
    rp_type_all{i} = rp_type;
    VMag_all{i} = v_mag;

    % Update the initial state for the next timestep
    state0 = states{end};
    s = EOSSeparator('pressure', 100000*Pascal, 'T', 273.15+22*Kelvin);
    sg = SeparatorGroup(s);
    sg.mode = 'moles';
    modelSparseAD.FacilityModel.SeparatorGroup = sg;

    % Display end of the current step
    disp(['Iteration ', num2str(i), ' out of ', num2str(num_iterations), ' completed']);
end


%% POST PROCESSING AND PLOTTING OF WELL DATA
Time_array = num2cell(repmat(86400, 1, num_iterations-1));
Time_array = cumsum(cell2mat(Time_array));
Time_array = [dt(1), Time_array]';  % dt(1) is the first element from the dt 
Time_array = num2cell(Time_array);

figure;
if ~isempty(own_wellSols{end})
    plotWellSols(own_wellSols, cell2mat(Time_array))
end


%% POST PROCESSING AND PLOTTING OF GRID DATA
% Adjust the saturations in own_states
for i = 1:numel(own_states)
    own_states{i}.s(abs(own_states{i}.s(:,1) - 1) < 1e-4, 1) = 0; % if oil saturation is 1 or close to 1, change it to 0
    own_states{i}.s(abs(own_states{i}.s(:,2)) < 1e-4, 2) = 1; % if gas saturation is 0 or close to 0, change it to 1
end
figure;
% Plot the results using plotToolbar
plotToolbar(modelSparseAD.G, own_states);

%% SAVE THE DATA PERMANENTLY FOR PLOTTING IN PYTHON
% Save WellSols to a MAT-file

baseDir = 'C:\Users\d.dorhjie\Documents\MATLAB\mrst-2023b\Bond Number\SPE\TAUS';

% save(fullfile(baseDir, 'WellSolsData_bond.mat'), 'own_wellSols');
% save(fullfile(baseDir, 'Schedule_bond.mat'), 'schedule');
% save(fullfile(baseDir, 'kor_all_bond.mat'), 'krr_all');
% save(fullfile(baseDir, 'own_states_data_bond.mat'), 'own_states');
% save(fullfile(baseDir, 'IFT_data_bond.mat'), 'IFT_C');
% save(fullfile(baseDir, 'dhro_bond.mat'), 'drho_all');
% save(fullfile(baseDir, 'Bond_number_bond.mat'), 'Bond_all');
% save(fullfile(baseDir, 'Capillary_number_bond.mat'), 'Capillary_all');
% save(fullfile(baseDir, 'RP_type_bond.mat'), 'rp_type_all');
% save(fullfile(baseDir, 'Vmag_bond.mat'), 'VMag_all');

save(fullfile(baseDir, 'WellSolsData_bond4.mat'), 'own_wellSols');
save(fullfile(baseDir, 'Schedule_bond4.mat'), 'schedule');
save(fullfile(baseDir, 'kor_all_bond4.mat'), 'krr_all');
save(fullfile(baseDir, 'own_states_data_bond4.mat'), 'own_states');
save(fullfile(baseDir, 'IFT_data_bond4.mat'), 'IFT_C');
save(fullfile(baseDir, 'dhro_bond4.mat'), 'drho_all');
save(fullfile(baseDir, 'Bond_number_bond4.mat'), 'Bond_all');
save(fullfile(baseDir, 'Capillary_number_bond4.mat'), 'Capillary_all');
save(fullfile(baseDir, 'RP_type_bond4.mat'), 'rp_type_all');
save(fullfile(baseDir, 'Vmag_bond4.mat'), 'VMag_all');


% End the timer and display the elapsed time
elapsedTime = toc/60;
fprintf('Total simulation time: %.2f minutes\n', elapsedTime);
% In this example:
%% COREY RELATIVE PERM FUNCTION WITHOUT EFFECT OF CAPILLARY NUMBER 
function fn = coreyPhaseRelpermAD_Des(n, sr, kwm, sr_tot)
    if nargin < 1
        n = 1;
    end
    if nargin < 2                                      
        sr = 0;
    end
    if nargin < 3
        kwm = 1;
    end
    if nargin < 4
        sr_tot = sr;
    end
    fn = @(s, varargin) coreyRelperm_Des(s, n, sr, kwm, sr_tot);  % commened by Desmond
end


% corey fucntion by the mrst team
function kr = coreyRelperm_Des(s, n, sr, kwm, sr_tot)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end

    kr = kwm*sat.^n;
end


%% COREY RELATIVE PERMEABILITY FUNCTON WITH EFFECT OF CAPILLARY NUMBER
function fn = coreyPhaseRelpermAD_DesNc(n, sr, kwm, sr_tot,Nc)  
    if nargin < 1
        n = 1;
    end
    if nargin < 2
        sr = 0;
    end
    if nargin < 3
        kwm = 1;
    end
    if nargin < 4
        sr_tot = sr;
    end
    fn = @(s, varargin) coreyRelperm_DesNc(s, n, sr, kwm, sr_tot,Nc); % inputed by Desmond
end

% corey function with dependance on capillary number as an exponent (this was just for try)
function kr = coreyRelperm_DesNc(s, n, sr, kwm, sr_tot, Nc)
    den = 1 - sr_tot;
    sat = ((s - sr) / den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end

    % if exist('Nc', 'var')
    %     kr = kwm * sat.^Nc;
    % else
    %     kr = kwm * sat.^n;
    % end
    % np = 1; % Assuming np is a constant
    % krbp = kwm * sat.^n;
    % krmp = sat;
    % kr = Nc.^(1/np) * krbp + (1 - Nc.^(1/np)) * krmp;
    %%
    mp = 1;
    % n1p = 1;
    % n2p = 1;
    % np = n1p * sat.^n2p;
    np = 1;
    Xp = 1 - exp(-mp .*Nc); %maybe put Tau here
    krbp = kwm * sat.^n;
    krmp = (sat - Xp*sr) ./ (1- Xp*sr);  % the miscible curve (normally X shaped RP)
    % krmp = sat;
    kr = Nc.^(1/np) * krbp + (1 - Nc.^(1/np)) * krmp;
end

%% COREY RELATIVE PERMEABILITY FUNCTON WITH EFFECT OF BOND NUMBER
function fn = coreyPhaseRelpermAD_DesNb(n, sr, kwm, sr_tot,Nb,Xp)  
    if nargin < 1
        n = 1;
    end
    if nargin < 2
        sr = 0;
    end
    if nargin < 3
        kwm = 1;
    end
    if nargin < 4
        sr_tot = sr;
    end
    fn = @(s, varargin) coreyRelperm_DesNb(s, n, sr, kwm, sr_tot,Nb,Xp); % inputed by Desmond
end

% % corey function with dependance on capillary number as an exponent (this was just for try)
% function kr = coreyRelperm_DesNb(s, n, sr, kwm, sr_tot, Nb,Xp)
%     den = 1 - sr_tot;
%     sat = ((s - sr) / den);
%     if isa(sat, 'ADI')
%         sat.val = max(min(sat.val, 1), 0);
%     else
%         sat = max(min(sat, 1), 0);
%     end
% 
%     % if exist('Nc', 'var')
%     %     kr = kwm * sat.^Nc;
%     % else
%     %     kr = kwm * sat.^n;
%     % end
%     % np = 1; % Assuming np is a constant
%     % krbp = kwm * sat.^n;
%     % krmp = sat;
%     % kr = Nc.^(1/np) * krbp + (1 - Nc.^(1/np)) * krmp;
%     %%
%     % mp = 1;
%     % n1p = 1;
%     % n2p = 1;
%     % np = n1p * sat.^n2p;
%     np = 2; %%
%     % Xp = 1 - exp(-mp .*Nb); %% %maybe put Tau here
%     krbp = kwm * sat.^n;
%     krmp = (sat - Xp*sr) ./ (1- Xp*sr);  % the miscible curve (normally X shaped RP)
%     % krmp = sat;
%     kr = Nb.^(1/np) * krbp + (1 - Nb.^(1/np)) * krmp;
% end
% corey function with dependance on capillary number as an exponent (this was just for try)
function kr = coreyRelperm_DesNb(s, n, sr, kwm, sr_tot, Nb,Xp)
    den = 1 - sr_tot;
    sat = ((s - sr) / den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
   
    % % coefficients
    % n1 =1;     %choosen for a purpose
    % n2 = 0.5;  %choosen for a purpose
    % np = n1*sat.^n2;
    np =1;

    % real relative permeability
    krbp = kwm * sat.^n;

    % x shaped relative permeability
    krmp = (sat - Xp*sr) ./ (1- Xp*sr);  % the miscible curve (normally X shaped RP)
    
    krmp = max(min(krmp, 1), 0);

    % outputed relative permeability
    % kr = Nb.^(1./np) .* krbp + (1 - Nb.^(1./np)) .* krmp;
    kr = Nb.^(1/np) * krbp + (1 - Nb.^(1/np)) * krmp;
end