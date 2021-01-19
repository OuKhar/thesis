% script to run macroeconomic DSGE models in dynare in oskari jakonen master's thesis
function []=run(model_selection)

if ~exist('model_selection','var')
    % if no input model is set, set default to the final model run
    model_selection = 'oskari';
end

clc;
close all;

% adjust path to where the dynare mod files are
cd([cd '/dynare_files']); 

if isequal(model_selection, 'mazelis')
% run dynare file for mazelis model with shadow banking
    dynare oskari_mazelis15.mod
    load('oskari_mazelis15_results.mat', 'M_', 'oo_')
elseif isequal(model_selection, 'gk')
% run dynare file for gertler-karadi model (without shadow banking)
    dynare oskari_gk11.mod
    load('oskari_gk11_results.mat', 'M_', 'oo_')
elseif isequal(model_selection, 'oskari')
% run dynare file for mazelis model with shadow banking altered to chinese
% economy
    dynare oskari_dsge_model.mod
    load('oskari_dsge_model_results.mat', 'M_', 'oo_') 
elseif isequal(model_selection, 'build')
% build the shadow banking model IRF's into single graphs 
    dynare oskari_dsge_model.mod
    load('oskari_dsge_model_results.mat', 'oo_')
    irf_oskari=oo_.irfs;
    save irf_oskari

    dynare oskari_mazelis15.mod
    load('oskari_mazelis15_results.mat', 'oo_')
    irf_maz=oo_.irfs;
    save irf_maz
    close all;

    load irf_oskari
    load irf_maz
    % define shocks
    ending_cell={'_e_A','_e_I'};

    for ii=1:length(ending_cell)
    % time horizon, max 40 periods defined in .mod files
    horizon=1:1:40;
    %define observed variables
    var={'y','c','int','k','inv','pi','sb','ssb','lab'};
    figure
    for jj=1:length(var)
        % 3x3 matrix to show the graphs for 9 vars
        subplot(3,3,jj)
        eval(['irf_oskari.' var{1,jj},ending_cell{1,ii}]);
        eval(['irf_maz.' var{1,jj},ending_cell{1,ii}]);
        hold on
        plot(horizon,[eval(['irf_oskari.' var{1,jj},ending_cell{1,ii}])],'--b',horizon, ...
            [eval(['irf_maz.' var{1,jj},ending_cell{1,ii}])],'-k','LineWidth',1);
    title([var{1,jj}])
    end
    end
else
    disp('no proper run argument specified')
    model_selected = ['run argument selected: ', model_selection];
    disp(model_selected)
    disp('run with "gk", "mazelis", "oskari" or "build" (check README)')
end

if isequal(model_selection,'mazelis') || isequal(model_selection,'gk') ...
        || isequal(model_selection,'oskari')
    % basic checks for the selected model
    disp('-----------------------')
    disp('performing basic checks')
    disp('-----------------------')
    disp('consumption of output')
    consumption =                                                           ...
    (exp(oo_.steady_state(strmatch('c',M_.endo_names))))                    ...
        /(exp(oo_.steady_state(strmatch('ym',M_.endo_names))));
    disp(consumption);
    disp('-----------------------')
    disp('labour of output')
    labour =                                                                ...
    (exp(oo_.steady_state(strmatch('lab',M_.endo_names))))                  ...
        /(exp(oo_.steady_state(strmatch('ym',M_.endo_names))));
    disp(labour);
end
    % shadow banking specific checks
if isequal(model_selection,'mazelis') || isequal(model_selection,'oskari')
    disp('-----------------------')
    disp('fund shares of savings')
    fundshares =                                                        ...
    (exp(oo_.steady_state(strmatch('fs',M_.endo_names))))               ...
        /(exp(oo_.steady_state(strmatch('b',M_.endo_names))));
    disp(fundshares);
    disp('-----------------------')
    disp('bank deposits of savings')
    consumption_from_output =                                           ...
    (exp(oo_.steady_state(strmatch('d ',M_.endo_names))))               ...
        /(exp(oo_.steady_state(strmatch('b',M_.endo_names))));
    disp(consumption_from_output);
end
% back to root
cd([cd '/..']);
end