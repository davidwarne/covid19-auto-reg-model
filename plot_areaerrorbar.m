% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one and two standard deviations;      %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval;              %
%                   if 'credint'    95% and 90% Bayesian credible interval%
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
%   Modified by: David Warne
%   E-mail:  david.warne@qut.edu.au
%   Data:    01/04/2020
% ----------------------------------------------------------------------- %
function [options] = plot_areaerrorbar(data,options)

    % Default options
    if(nargin<2)
        options.handle     = figure(1);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        options.alpha      = 0.5;
        options.line_width = 1;
        options.error      = 'std';
    end
    
    if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
    options.x_axis = options.x_axis(:);
    
    % Computing the mean and standard deviation of the data matrix
    data_mean = mean(data,1);
    data_std  = std(data,0,1);
    
    % Type of error plot
    switch(options.error)
        case 'std', error = data_std;
        case 'sem', error = (data_std./sqrt(size(data,1)));
        case 'var', error = (data_std.^2);
        case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
    
    end

    if strcmp(options.error,'credint')
        % Plotting the result
        figure(options.handle);
        % computing quantiles for 99% 95% 90% credible intervals
        Q = quantile(data,[0.5;0.25;0.75;0.05;0.95])

        hold on; % added DJW
        x_vector = [options.x_axis', fliplr(options.x_axis')];
        patch95 = fill(x_vector, [Q(5,:),fliplr(Q(4,:))], options.color_area); 
        set(patch95, 'edgecolor', 'none'); 
        set(patch95, 'FaceAlpha', options.alpha);
        patch90 = fill(x_vector, [Q(3,:),fliplr(Q(2,:))], options.color_line);
        set(patch90, 'edgecolor', 'none');
        set(patch90, 'FaceAlpha', options.alpha*1.25);
        %hold on;
        %plot(options.x_axis, Q(1,:), 'color', options.color_line, ...
        %    'LineWidth', options.line_width);
        hold off;
    else        
        % Plotting the result
        figure(options.handle);
        hold on; % added DJW
        x_vector = [options.x_axis', fliplr(options.x_axis')];
        patch = fill(x_vector, [data_mean+2*error,fliplr(data_mean-2*error)], options.color_area); 
        set(patch, 'edgecolor', 'none'); 
        set(patch, 'FaceAlpha', options.alpha);
        patch2 = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
        set(patch2, 'edgecolor', 'none');
        set(patch2, 'FaceAlpha', options.alpha);
        hold on;
        plot(options.x_axis, data_mean, 'color', options.color_line, ...
            'LineWidth', options.line_width);
        hold off;
    end
end
