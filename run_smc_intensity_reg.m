%% Bayesian analysis of COVID-19 pandemic
%  Utilises data provided by from various sources (including Johns-Hopkins University
%   DXY, National Health Commisson (China), Protezione Civile (Italy), 
%   and https://www.corona-data.ch/ (Switzerland)).
%  
% Simulation based model extending a stochastic SIR model with latent infectous
% populatio and regulatory mechnisms.
%
%  Inference of model parameters is perfromed using Approximate Bayesian Computation
%  Within the Sequential Monte Carlo framework of Drovandi and Pettit (2011).
%  
% Authors:
%     Christopher Drovandi (c.drovandi@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%
%     David J. Warne (david.warne@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Algorithm Initialisation
%
% Initialise Random number generator for reproducibility
%clear global
% country_id defined within the current workspace to select the region to perform
% inference on.
rng(1337,'twister')

%%
region_id = 'Italy'; % 
iso_id = 'ITA';
province_id = 'total';
results.name = sprintf('%s%s',region_id,province_id);
results.ISO3166alpha3 = iso_id;

%% data store
DATA_DIR = '../'

% load transition matrix for confirmed cases by country
opts = detectImportOptions([DATA_DIR,'COVID19data/data-raw/covid19_sorted.csv']);
for i=1:length(opts.VariableNames)
    switch opts.VariableNames{i}
        case 'alpha3'
            opts.VariableTypes{i} = 'categorical';
        case 'Country_Regions'
            opts.VariableTypes{i} = 'categorical';
        case 'Province_State'
            opts.VariableTypes{i} = 'categorical';
        case 'date'
            opts.VariableTypes{i} = 'datetime';
        otherwise
            opts.VariableTypes{i} = 'double';
    end
end
T_d = readtable([DATA_DIR,'COVID19data/data-raw/covid19_sorted.csv'],opts);

% load population table
opts = detectImportOptions([DATA_DIR,'COVID-19_ISO-3166/full_list.csv']);
for i=1:length(opts.VariableNames)
    switch opts.VariableNames{i}
        case 'alpha3'
            opts.VariableTypes{i} = 'categorical';
        case 'Province_State'
            opts.VariableTypes{i} = 'categorical';
        case 'population'
            opts.VariableTypes{i} = 'double';
    end
end
P_d = readtable([DATA_DIR,'COVID-19_ISO-3166/full_list.csv'],opts);

%% extract region of interest
I = T_d.alpha3 == iso_id & T_d.Province_State == province_id;
Data.D = T_d.deaths(I);
Data.R = T_d.recovered(I);
Data.C = T_d.active(I);
J = P_d.alpha3 == iso_id & P_d.Province_State == province_id;
Data.P = P_d.population(J);

% check simulation is viable
if Data.P <= 0 || sum(Data.C > 0) == 0 || P_d.alpha3(J) == 'cruise'
    % unpopulated region or no cases of COVID-19 reported
    return
end

%% Synthetic data for identifiability test
%Data = simuldata_reg(Data,[0.01,41,0.1,0.05,0.1,1,1/2,1])

%% Back testing (only fit up to day T-test_days and use remainer to validate)
%test_data = Data;
test_data = Data;
test_days = 5;
T = length(Data.C);
Data.C(T-test_days:T) = [];
Data.R(T-test_days:T) = [];
Data.D(T-test_days:T) = [];

%% set-up ABC-SMC 
% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.75; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.005;

% define prior
prior.num_params = 8; % [\alpha_0,\alpha,\beta,\gamma,\alpha_tau]
prior.p1 = [0.0  ,0.0, 0.0,0.0,0.0,0.0,0.0,0];
prior.p2 = [2.0 ,100.0, 1.0, 1.0, 1.0, 1.0,2,2];
prior.sampler = @() [unifrnd(prior.p1,prior.p2)]; 
prior.pdf = @(theta) prod([unifpdf(theta,prior.p1,prior.p2)]);
prior.trans_f = @(theta) [theta];
prior.trans_finv = @(theta) [theta];
% user functions

% stochastic simulation
sim_func = @simuldata_reg;
%discrepancy metric
dist_func = @(S_d,S_s) sqrt(sum((S_d(:) - S_s(:)).^2));
% summary statistic
smry_func = @smry;

%% run SMC sampler
[part_vals, part_sim, part_s, sims,epsilon_t,p_acc_t] = smc_abc_rw_par(Data,...
                                       sim_func,dist_func,smry_func,prior,N,...
                                       epsilon_final,a,c,p_acc_min);
%% save results
results.part_vals = part_vals;
results.part_sim = part_sim;
results.part_s = part_s;
results.sims = sims;
results.epsilon = epsilon_t;
results.p_acc_t = p_acc_t;
results.data = Data;
% save output
save(['results_smc_intensity_poisson_reg',results.name,'_',results.ISO3166alpha3,'.mat'],'results');
%% load output 
load(['results_smc_intensity_poisson_reg',results.name,'_',results.ISO3166alpha3,'.mat'],'results');

%% plot marginal posterior densities
figure(1);
lab = {'\alpha_0','\alpha','\beta','\gamma','\delta','\eta','n','\kappa'};
for i = 1:8
    subplot(2,4,i);
    ksdensity(results.part_vals(:,i));
    xlim([prior.p1(i),prior.p2(i)])
    xlabel(lab{i})
end
title(sprintf('%s %s',results.name,results.ISO3166alpha3))

%% plot samples of posterior predictve distribution against data
%figure;
T = length(results.data.C);
optsf = [];
optsf.handle = figure(3);
optsf.line_width = 1;
optsf.alpha = 0.5;
optsf.error = 'std';
optsf.color_area = [128 193 219]./255;    % Blue theme
optsf.color_line = [ 52 148 186]./255;
optsf = plot_areaerrorbar(results.part_sim(:,1:T),optsf)
optsf.color_area = [237,177,32]./255;    % yellow theme
optsf.color_line = [237,177,32]./(2*255);
optsf = plot_areaerrorbar(results.part_sim(:,T+1:2*T),optsf)
optsf.color_area = [243 169 114]./255;    % Orange theme
optsf.color_line = [236 112  22]./255;
optsf = plot_areaerrorbar(results.part_sim(:,2*T+1:3*T),optsf)
hold on
title(sprintf('Fit %s %s',results.name,results.ISO3166alpha3))
plot(results.data.C,'+k');
plot(results.data.R,'+k');
plot(results.data.D,'+k');
xlabel('time');
ylabel('counts')

pred_C = mean(results.part_sim(:,1:T))';
pred_R = mean(results.part_sim(:,T+1:2*T))';
pred_D = mean(results.part_sim(:,2*T+1:3*T))';

chisq_C = ((pred_C - results.data.C).^2)./var(results.part_sim(:,1:T))';
chisq_R = ((pred_R - results.data.R).^2)./var(results.part_sim(:,T+1:2*T))';
chisq_D = ((pred_D - results.data.D).^2)./var(results.part_sim(:,2*T+1:3*T))';
%
fit_T = table(pred_C,pred_R,pred_D,results.data.C,results.data.R,results.data.D,chisq_C,chisq_R,chisq_D,...
    'VariableNames',{'pred_C','pred_R','pred_D','obs_C','obs_R','obs_D','chisq_C','chisq_R','chisq_D'});

%% forward predictions for unobserved future
%figure;
pred_days = 10;
data_pred = results.data;
%T = length(test_data.C);
optsf = [];
optsf.handle = figure(4);
optsf.line_width = 1;
optsf.alpha = 0.5;
optsf.error = 'std';
% append 10 extra days for simulation
data_pred.C = [data_pred.C;zeros(pred_days,1)];
data_pred.D = [data_pred.D;zeros(pred_days,1)];
data_pred.R = [data_pred.R;zeros(pred_days,1)];
T = length(data_pred.C);
predsims = zeros(N,3*T);
for i=1:N
    D_s = sim_func(data_pred,results.part_vals(i,:));
    predsims(i,:) = smry(D_s);
end
optsf.color_area = [128 193 219]./255;    % Blue theme
optsf.color_line = [ 52 148 186]./255;
optsf = plot_areaerrorbar(predsims(:,1:T),optsf)
optsf.color_area = [237,177,32]./255;    % yellow theme
optsf.color_line = [237,177,32]./(2*255);
optsf = plot_areaerrorbar(predsims(:,T+1:2*T),optsf)
optsf.color_area = [243 169 114]./255;    % Orange theme
optsf.color_line = [236 112  22]./255;
optsf = plot_areaerrorbar(predsims(:,2*T+1:3*T),optsf)
hold on
title(sprintf('Predictions %s %s',results.name,results.ISO3166alpha3))
plot(results.data.C,'+k');
plot(results.data.R,'+k');
plot(results.data.D,'+k');
xlabel('time');
ylabel('counts')
save(['smc_intensity_poisson_reg_predictions',results.name,'_',results.ISO3166alpha3,'.fig'])
saveas(gcf,['smc_intensity_poisson_reg_predictions',results.name,'_',results.ISO3166alpha3,'.pdf'])

writetable(fit_T, ['smc_intensity_poisson_reg_predictions',results.name,'_',results.ISO3166alpha3,'.csv'])

%% Back testing validation plot forward predictions data
T = length(test_data.C);
optsf = [];
optsf.handle = figure(5);
optsf.line_width = 1;
optsf.alpha = 0.5;
optsf.error = 'std';
predsims = zeros(N,3*T);
for i=1:1000
   D_s = sim_func(test_data,results.part_vals(i,:));
   predsims(i,:) = smry(D_s);
end
optsf.color_area = [128 193 219]./255;    % Blue theme
optsf.color_line = [ 52 148 186]./255;
optsf = plot_areaerrorbar(predsims(:,1:T),optsf)
optsf.color_area = [237,177,32]./255;    % yellow theme
optsf.color_line = [237,177,32]./(2*255);
optsf = plot_areaerrorbar(predsims(:,T+1:2*T),optsf)
optsf.color_area = [243 169 114]./255;    % Orange theme
optsf.color_line = [236 112  22]./255;
optsf = plot_areaerrorbar(predsims(:,2*T+1:3*T),optsf)
hold on
title(sprintf('Predictions %s %s',results.name,results.ISO3166alpha3))
plot(results.data.C,'+k');
plot(results.data.R,'+k');
plot(results.data.D,'+k');
plot(test_data.C,'sk');
plot(test_data.R,'sk');
plot(test_data.D,'sk');
xlabel('time');
ylabel('counts')

pred_C = mean(predsims(:,1:T))';
pred_R = mean(predsims(:,T+1:2*T))';
pred_D = mean(predsims(:,2*T+1:3*T))';

chisq_C =  ((pred_C - test_data.C).^2)./var(predsims(:,1:T))';
chisq_R = ((pred_R - test_data.R).^2)./var(predsims(:,T+1:2*T))';
chisq_D = ((pred_D - test_data.D).^2)./var(predsims(:,2*T+1:3*T))';
%
pred_T = table(pred_C,pred_R,pred_D,test_data.C,test_data.R,test_data.D,chisq_C,chisq_R,chisq_D,...
   'VariableNames',{'pred_C','pred_R','pred_D','obs_C','obs_R','obs_D','chisq_C','chisq_R','chisq_D'});

save(['smc_intensity_poisson_reg_backtesting',results.name,'_',results.ISO3166alpha3,'.fig'])
saveas(gcf,['smc_intensity_poisson_reg_backtesting',results.name,'_',results.ISO3166alpha3,'.pdf'])

writetable(pred_T, ['smc_intensity_poisson_reg_backtesting',results.name,'_',results.ISO3166alpha3,'.csv'])
