%% Bayesian analysis of COVID-19 pandemic
%  Utilises data provided by from various sources (including Johns-Hopkins University
%   DXY, National Health Commisson (China), Protezione Civile (Italy), 
%   and https://www.corona-data.ch/ (Switzerland)).
%  
% post processing analysis and plotting
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: plotting functions desiged for processing all results from a run on the cluster.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

data_date_pred = '25jun';

data_date = '10jun';
pred_horizon = 0;
%model = 0;
model = 0;
% load population table

%% data store
DATA_DIR = './data/';
RES_DIR = ['./results_',data_date,'/T',num2str(pred_horizon),'/M',num2str(model)];

% load transition matrix for confirmed cases by country
opts = detectImportOptions([DATA_DIR,'covid19_sorted_',data_date,'.csv']);
% note: this is a bit verbose, but it is the most effective way I can get
% to ensure the order and number of columns is arbitrary
for i=1:length(opts.VariableNames)
    switch opts.VariableNames{i}
        case 'iso3'
            opts.VariableTypes{i} = 'categorical';
        case 'alpha3'
            opts.VariableTypes{i} = 'categorical';
        case 'Country_Region'
            opts.VariableTypes{i} = 'categorical';
        case 'Province_State'
            opts.VariableTypes{i} = 'categorical';
        case 'date'
            opts.VariableTypes{i} = 'datetime';
        otherwise
            opts.VariableTypes{i} = 'double';
    end
end
T_d = readtable([DATA_DIR,'covid19_sorted_',data_date,'.csv'],opts);
opts = detectImportOptions([DATA_DIR,'covid19_sorted_',data_date_pred,'.csv']);
for i=1:length(opts.VariableNames)
    switch opts.VariableNames{i}
        case 'iso3'
            opts.VariableTypes{i} = 'categorical';
        case 'alpha3'
            opts.VariableTypes{i} = 'categorical';
        case 'Country_Region'
            opts.VariableTypes{i} = 'categorical';
        case 'Province_State'
            opts.VariableTypes{i} = 'categorical';
        case 'date'
            opts.VariableTypes{i} = 'datetime';
        otherwise
            opts.VariableTypes{i} = 'double';
    end
end
T_d_pred = readtable([DATA_DIR,'covid19_sorted_',data_date_pred,'.csv'],opts);

% stochastic simulation (use @simuldata_reg_fD for death tally only version)
sim_func = @simuldata_reg_fA;
if model == 0
    sim_func = @simuldata_reg_fA;
elseif model == 1
    sim_func = @simuldata_reg_fC;
elseif model == 2
    sim_func = @simuldata_reg_fD;
elseif model == 3
    sim_func = @simuldata_reg;
elseif model == 4
    sim_func = @simuldata_reg_U;
end

%discrepancy metric
dist_func = @(S_d,S_s) sqrt(sum((S_d(:) - S_s(:)).^2));
% summary statistic
smry_func = @smry;

%% load output 
lab = {'\alpha_0','\alpha','\beta','\gamma','\delta','\eta','n','\kappa','-\log w_A'};
ii = 1;
for j=1:252
    try
        fprintf([RES_DIR,'/results_smc_intnsity_poisson_reg',num2str(j),'.mat\n']);
        load([RES_DIR,'/results_smc_intensity_poisson_reg',num2str(j),'.mat'],'results');
        results_array(ii) = results;
        ii = ii + 1;
    end
end

iso_id = []
for i=1:length(results_array)
    iso_id = [iso_id;results_array(i).ISO3166alpha3]
end

%% box plots
lab = {'\alpha_0','\alpha','\beta','\gamma','\delta','\eta','n','\kappa','-\log w_A'};
panel =1
for j=1:length(results_array(1).part_vals(1,1:9))
    X = [];
    G = [];
    for i=1:length(results_array)
        if results_array(i).ISO3166alpha3 ~= 'SMR'  && results_array(i).ISO3166alpha3 ~= 'BLR' && results_array(i).ISO3166alpha3 ~= 'GEO'
            X = [X;results_array(i).part_vals(1:500,j)]; 
            G = [G;repmat(results_array(i).ISO3166alpha3,length(results_array(i).part_s(1:500)),1)];
        end
    end
    figure(panel)
    if mod(j,3) == 0
        panel = panel +1;
    end
    subplot(3,1,mod(j+2,3)+1)
    boxplot(X,G,'PlotStyle','compact','symbol','');
    xlabel('Country Code (ISO-3166 alpha-3)');
    ylabel(['$',lab{j}, ' \sim p(',lab{j},' | \mathcal{D})$'])
    ytickformat('%0.2f')
end

%% point estimates
% compute lists for sorting and ranking
C_list = zeros(length(results_array),1);
R_list = zeros(length(results_array),1);
D_list = zeros(length(results_array),1);

max_dC_list = zeros(length(results_array),1);
max_dR_list = zeros(length(results_array),1);
max_dD_list = zeros(length(results_array),1);
dC_list = zeros(length(results_array),1)
for i=1:length(results_array)
    C_list(i) = results_array(i).data.C(end) + results_array(i).data.R(end) + results_array(i).data.D(end);
    R_list(i) = results_array(i).data.R(end);
    D_list(i) = results_array(i).data.D(end);
    max_dC_list(i) = max(diff(results_array(i).data.C+ results_array(i).data.R+results_array(i).data.D));
    max_dR_list(i) = max(diff(results_array(i).data.R));
    max_dD_list(i) = max(diff(results_array(i).data.D));
    dC_list(i) = diff(results_array(i).data.C(end-1:end)+ results_array(i).data.R(end-1:end)+results_array(i).data.D(end-1:end));
end

point_ests = zeros(length(results_array),length(results_array(1).part_vals(1,:)));
CI.low = zeros(length(results_array),length(results_array(1).part_vals(1,:)));
CI.high = zeros(length(results_array),length(results_array(1).part_vals(1,:)));
CI.c95 = zeros(length(results_array),length(results_array(1).part_vals(1,:)));

% compute estimates
for i=1:length(results_array)
    point_ests(i,1:length(results_array(i).point_est)) = results_array(i).point_est;
    [~,I] = sort(results_array(i).part_s);
    temp = results_array(i).part_vals(I,:);
    temp(:,9) = 10.^(-temp(:,9));
    CI.low(i,:) = quantile(temp(1:500,:),0.025);
    CI.high(i,:) = quantile(temp(1:500,:),0.975);
    CI.c95(i,:) = std(results_array(i).part_vals)*1.96/sqrt(size(results_array(i).part_vals,1));
end

output_table = table(iso_id,...
                     point_ests(:,1),CI.low(:,1),CI.high(:,1),CI.c95(:,1),...
                     point_ests(:,2),CI.low(:,2),CI.high(:,2),CI.c95(:,2),...
                     point_ests(:,3),CI.low(:,3),CI.high(:,3),CI.c95(:,3),...
                     point_ests(:,4),CI.low(:,4),CI.high(:,4),CI.c95(:,4),...
                     point_ests(:,5),CI.low(:,5),CI.high(:,5),CI.c95(:,5),...
                     point_ests(:,6),CI.low(:,6),CI.high(:,6),CI.c95(:,6),...
                     point_ests(:,7),CI.low(:,7),CI.high(:,7),CI.c95(:,7),...
                     point_ests(:,8),CI.low(:,8),CI.high(:,8),CI.c95(:,8),...
                     point_ests(:,9),CI.low(:,9),CI.high(:,9),CI.c95(:,9),...
                     C_list,R_list,D_list,... 
 'VariableNames',{'iso_id',...
                    'alpha_0','alpha_0_CrIL','alpha_0_CrIU','alpha_0_CI',...
                    'alpha',  'alpha_CrIL','alpha_CrIU','alpha_CI',...
                    'beta','beta_CrIL','beta_CrIU','beta_CI',...
                    'gamma','gamma_CrIL','gamma_CrIU','gamma_CI',...
                    'delta','delta_CrIL','delta_CrIU','delta_CI',...
                    'eta','eta_CrIL','eta_CrIU','eta_CI',...
                    'n','n_CrIL','n_CrIU','n_CI',...
                    'kappa','kappa_CrIL','kappa_CrIU','kappa_CI',...
                    'wA','wA_CrIL','wA_CrIU','wA_CI',...
                    'C_T','R_T','D_T'});
results_table = table(point_ests(:,1),point_ests(:,2),point_ests(:,3),...
                       point_ests(:,4),point_ests(:,5),point_ests(:,6),...
                       point_ests(:,7),point_ests(:,8),10.^(-point_ests(:,9)),10.^(-point_ests(:,10)),10.^(-point_ests(:,11)),...
                       C_list,R_list,D_list,max_dC_list,max_dR_list,max_dD_list,iso_id,...
 'VariableNames',{'alpha_0','alpha','beta','gamma','delta','eta','n','kappa','wA','wD','wR','C_T','R_T','D_T','dC','dR','dD','iso_id'});

fid = 1
%% plot scatter plots of point estimates
rlist = {'USA','BRA','RUS','GBR','IND','ESP','PER','FRA'}
glist = {'CHN','KOR','AUS','NZL','ITA','DEU'}
figure;
subplot(1,6,1)
semilogy(results_table.n,results_table.log_wA,'+','Color',[1,1,1]/2)
hold on
semilogy(results_table.n(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),'+g')
semilogy(results_table.n(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),'og')
semilogy(results_table.n(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),'+r')
semilogy(results_table.n(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),'or')
text(results_table.n(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),results_table.iso_id(ismember(results_table.iso_id, glist)));
text(results_table.n(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),results_table.iso_id(ismember(results_table.iso_id, rlist)));
xlabel('$n$');ylabel('$w_A$');

subplot(1,6,2)
semilogy(results_table.gamma,results_table.log_wA,'+','Color',[1,1,1]/2)
hold on
semilogy(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),'+g')
semilogy(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),'og')
semilogy(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),'+r')
semilogy(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),'or')
text(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),results_table.iso_id(ismember(results_table.iso_id, glist)));
text(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),results_table.iso_id(ismember(results_table.iso_id, rlist)));
xlabel('$\gamma$');ylabel('$w_A$');

subplot(1,6,3)
loglog(results_table.kappa,results_table.log_wA,'+','Color',[1,1,1]/2)
hold on
loglog(results_table.kappa(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),'+g')
loglog(results_table.kappa(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),'og')
loglog(results_table.kappa(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),'+r')
loglog(results_table.kappa(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),'or')
text(results_table.kappa(ismember(results_table.iso_id, glist)),results_table.log_wA(ismember(results_table.iso_id, glist)),results_table.iso_id(ismember(results_table.iso_id, glist)));
text(results_table.kappa(ismember(results_table.iso_id, rlist)),results_table.log_wA(ismember(results_table.iso_id, rlist)),results_table.iso_id(ismember(results_table.iso_id, rlist)));
xlabel('$\kappa$');ylabel('$w_A$');

subplot(1,6,4)
semilogy(results_table.n,results_table.kappa,'+','Color',[1,1,1]/2)
hold on
semilogy(results_table.n(ismember(results_table.iso_id, glist)),results_table.kappa(ismember(results_table.iso_id, glist)),'+g')
semilogy(results_table.n(ismember(results_table.iso_id, glist)),results_table.kappa(ismember(results_table.iso_id, glist)),'og')
semilogy(results_table.n(ismember(results_table.iso_id, rlist)),results_table.kappa(ismember(results_table.iso_id, rlist)),'+r')
semilogy(results_table.n(ismember(results_table.iso_id, rlist)),results_table.kappa(ismember(results_table.iso_id, rlist)),'or')
text(results_table.n(ismember(results_table.iso_id, glist)),results_table.kappa(ismember(results_table.iso_id, glist)),results_table.iso_id(ismember(results_table.iso_id, glist)));
text(results_table.n(ismember(results_table.iso_id, rlist)),results_table.kappa(ismember(results_table.iso_id, rlist)),results_table.iso_id(ismember(results_table.iso_id, rlist)));
xlabel('$n$');ylabel('$\kappa$');

subplot(1,6,5)
semilogy(results_table.gamma,results_table.kappa,'+','Color',[1,1,1]/2)
hold on
semilogy(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.kappa(ismember(results_table.iso_id, glist)),'+g')
semilogy(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.kappa(ismember(results_table.iso_id, glist)),'og')
semilogy(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.kappa(ismember(results_table.iso_id, rlist)),'+r')
semilogy(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.kappa(ismember(results_table.iso_id, rlist)),'or')
text(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.kappa(ismember(results_table.iso_id, glist)),results_table.iso_id(ismember(results_table.iso_id, glist)));
text(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.kappa(ismember(results_table.iso_id, rlist)),results_table.iso_id(ismember(results_table.iso_id, rlist)));
xlabel('$\gamma$');ylabel('$\kappa$');

subplot(1,6,6)
plot(results_table.gamma,results_table.n,'+','Color',[1,1,1]/2)
hold on
plot(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.n(ismember(results_table.iso_id, glist)),'+g')
plot(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.n(ismember(results_table.iso_id, glist)),'og')
plot(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.n(ismember(results_table.iso_id, rlist)),'+r')
plot(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.n(ismember(results_table.iso_id, rlist)),'or')
text(results_table.gamma(ismember(results_table.iso_id, glist)),results_table.n(ismember(results_table.iso_id, glist)),results_table.iso_id(ismember(results_table.iso_id, glist)));
text(results_table.gamma(ismember(results_table.iso_id, rlist)),results_table.n(ismember(results_table.iso_id, rlist)),results_table.iso_id(ismember(results_table.iso_id, rlist)));
xlabel('$\gamma$');ylabel('$n$');


%% posterior predictive sampling
J = find(ismember(iso_id,{'CHN','KOR','TWN','SGP','AUS','NZL','USA','ESP','DEU','ITA','CHE','IRN','FRA','GBR','TUR'}));
fid = 1;
ii =1;
h = figure;
for j=J'
    results = results_array(j);
    j
    
    %% resampling and forward predictions for unobserved future
    pred_days = 0;
    data_pred = results.data;
    optsf = [];
    optsf.handle = figure(fid+1);
    optsf.line_width = 1;
    optsf.alpha = 0.5;
    optsf.error = 'credint';
    % append 10 extra days for simulation
    data_pred.C = [data_pred.C;zeros(pred_days,1)];
    data_pred.D = [data_pred.D;zeros(pred_days,1)];
    data_pred.R = [data_pred.R;zeros(pred_days,1)];
    T = length(data_pred.C);
    N = size(results.part_vals,1);
    predsims = zeros(N,3*T);
    for i=1:N
        %D_s = sim_func(results.data,results.point_est);
        D_s = sim_func(data_pred,results.part_vals(i,:));
        predsim = smry(D_s);
        predsims(i,:) = predsim;
    end
    
    time_seq = T_d.date(T_d.alpha3 == iso_id(j) & T_d.Province_State == 'total')';
    col = [253,191,111; 178,223,138;251,154,153] ./255
    figure(fid+1);
    subplot(1,3,1)
    pred_dC = diff(T_d_pred.confirmed(T_d_pred.alpha3 == results.ISO3166alpha3 & T_d_pred.Province_State == 'total'));
    pred_dR = diff(T_d_pred.recovered(T_d_pred.alpha3 == results.ISO3166alpha3 & T_d_pred.Province_State == 'total'));
    pred_dD = diff(T_d_pred.deaths(T_d_pred.alpha3 == results.ISO3166alpha3 & T_d_pred.Province_State == 'total'));
    bar(time_seq(2:end),[diff(results.data.D + results.data.R + results.data.C)],'FaceColor',col(1,:));
    I = find(T_d.alpha3 == results.ISO3166alpha3);
    title(sprintf('%s (%s)',T_d.Country_Region(I(1)),results.ISO3166alpha3));
    
     hold on
     bar([time_seq(end)+1:time_seq(end)+pred_days],pred_dC(end-pred_days+1:end),'FaceColor','r');
     ylabel('Daily number of cases');
     
     subplot(1,3,2)
     bar(time_seq(2:end),[diff(results.data.R)],'FaceColor',col(2,:));
     hold on
     ylabel('Daily number of recoveries');
     
     subplot(1,3,3)
     bar(time_seq(2:end),[diff(results.data.D)],'FaceColor',col(3,:));
     hold on
     ylabel('Daily number of deaths');
      
    optsf.x_axis = [time_seq(2:end), time_seq(end)+1:time_seq(end)+pred_days];
    optsf.color_area = [166,206,227]./255;    % Blue theme
    optsf.color_line = [ 31,120,180]./255;
    subplot(1,3,1)
    optsf = plot_areaerrorbar([diff(( predsims(:,1:T) + predsims(:,T+1:2*T) + predsims(:,2*T+1:3*T))')'],optsf)   

    optsf.color_area = [166,206,227]./255;    % Blue theme
    optsf.color_line = [ 31,120,180]./255;
    subplot(1,3,2)
    optsf = plot_areaerrorbar([diff(( predsims(:,T+1:2*T))')'],optsf)   
    
    optsf.color_area = [166,206,227]./255;    % Blue theme
    optsf.color_line = [ 31,120,180]./255;
    subplot(1,3,3)
    optsf = plot_areaerrorbar([diff(( predsims(:,2*T+1:3*T))')'],optsf) 

    ii = ii + 1;
    t = find(results.data.C + results.data.R + results.data.D >= 100);
    if length(t) > 2
        subplot(1,3,1);
        xlim([time_seq(t(2)), time_seq(end)+pred_days]);
        xtickformat('MMM-dd')
        subplot(1,3,2);
        xlim([time_seq(t(2)), time_seq(end)+pred_days]);
        xtickformat('MMM-dd')
        subplot(1,3,3);
        xlim([time_seq(t(2)), time_seq(end)+pred_days]);
        xtickformat('MMM-dd')
    end
    box on;
    set(gcf,'Units','centimeters','OuterPosition', [15.0707,10.1812,20.3835,8.3608])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = [20,6];
    
    fid = fid + 1
end

figure;
h=heatmap(corr(results_table{:,1:12},'Type','Spearman'));

h.YDisplayLabels = repmat({''}, size(h.YData));
h.XDisplayLabels = repmat({''}, size(h.XData));

