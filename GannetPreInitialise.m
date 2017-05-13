function MRS_struct=GannetPreInitialise(MRS_struct)

% Some of these parameters will be parsed from the data headers

% Acquisition Parameters
    MRS_struct.p.sw = []; % parsed from header
    MRS_struct.p.npoints = []; % parsed from header
    MRS_struct.p.TR = []; % parsed from header
    MRS_struct.p.TE = []; % parsed from header
    MRS_struct.p.LarmorFreq = []; % parsed from header
    MRS_struct.p.Nwateravg = []; % parsed from header
    MRS_struct.p.target = 'GABAGlx'; % options are 'GABA', 'Glx', 'GABAGlx' and 'GSH'
    MRS_struct.p.ONOFForder = 'offfirst'; % options are 'onfirst' or 'offfirst'
    MRS_struct.p.Water_Positive = 1; % for Philips MOIST ws, set to 0
    % Siemens header information differs between versions; switch for different versions
    MRS_struct.p.Siemens_type = 1; % 1 = TIM TRIO WIP; 2 = Near seq; 3 = Skyra WIP; 4 = Prisma (VD13C); 5 = Prisma (Minnesota); 6 = Jamie's VE11B (Jena)
    
% Analysis Parameters
    MRS_struct.p.LB = 3;
    %MRS_struct.p.ZeroFillTo = []; % zero-fill to obtain nominal spectral resolution of 0.061 Hz/point
    MRS_struct.p.AlignTo = 'SpecReg'; % options are 'SpecReg' (default and recommended), 'Cr', 'Cho', 'NAA', 'H20', 'CrOFF'
    
% Output Parameters
    MRS_struct.p.mat = 0;  % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.sdat = 0; % 1 = YES, save MRS_struct as .sdat file
    
end
