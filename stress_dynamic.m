% Copyright David F. Gold (dfg42@cornell.edu)
% Modified by Lillian B.J. Lau (lbl59@cornell.edu) 
% Scenario-Tunable Regional Synthetic Streamflow (STRESS) Generator
%
% Q_historical is the historical streamflow cell array
% num_syn_years is the number of synthetic years to generate
% p is % of lowest streamflows that are selected
%   Change the magnitude of p to determine if you want more extreme
%   floods or more extreme droughts.
%   If p = 0, outputs stationary synthetic flows.
% n is a vector containing the number of copies of the lowest p% of flows
% m is a vector containing the number of copies of the highest (100-p)% 
%   of flows
%
% The length of n and m are the number of perturbation periods to be 
% generated throughout the synthetic record
%
% Each period has a different flood/drought risk 
%

function Qs = stress_dynamic(Q_historical, num_syn_years, p, n, m)
    % Input error checking
    if iscell(Q_historical)
        npoints = length(Q_historical);
        % npoints is equal to the number of streamflow stations to generate
        % synthetic flows for
    else
        error(['Q_historical must be a cell array ' ...
               'containing one or more 2-D matrices.']);
    end
    
    nQ_historical = length(Q_historical{1}(:,1));
    % nQ_historical is the number of years in the historical record
    % at each streamflow station
    
    for i=2:npoints
        if length(Q_historical{i}(:,1)) ~= nQ_historical
            error('All matrices in Q_historical must be the same size.');
        end
    end
    % Input error checking ends
    num_syn_years = num_syn_years+1; % adjust for the new corr technique
    weekly_mean = zeros(52,1);
    weekly_stdev = zeros(52,1);
    period_length = 1;  % default value
    Random_Matrix = zeros(num_syn_years,52);    % matrix M
    % Random Matrix will keep the random number
    % which will be the (historical) year from which we will bootstrap
    % values from the historical record at each streamflow site
    % for each week in the synthetic record (its the random
    % component to this process) note: this is done at
    % the beginning so that the random number is the same
    % for each streamflow site in a given week
    
    if p > 0      
        period_length = num_syn_years/length(n);
    end
    % period length is the number of years each
    % perturbation period lasts (i.e. - if n & m have five values
    % then there are five separate perturbation periods and
    % each period lasts 1/5 of the synthetic record
    
    for yr = 1:num_syn_years 
        period = ceil(yr/period_length);
        % determine which perturbation period we are in, this affects the 
        % length of the record (see immediately below)
        if p == 0
        	nQ = nQ_historical;
            % two input arguments give you the synthetic flows 
            % from a stationary record
    	
        elseif p > 0
            nQ = nQ_historical + ceil((p*nQ_historical))*n(period) + (nQ_historical - ceil((p*nQ_historical)))*m(period);
            % runs in perturbed mode - the perturbations actually increase
            % the length of the record so we need to adjust the nQ value
            % for bootstrapping - this is the length of the 'perturbed
            % historical' record. This adds n copies of the lowest p% of
            % streamflow years to the length of the historical record (so 10% of 81
            % years is 9 years (round up), and 5 times more likely would add 40 years)
            % and m copies of the highest (1 - p)% of streamflow years, (so
            % 90% of 81 years is 73 years, times 2 would be 146 years,
            % giving us a new record length of 81 + 40 + 146 or 267 years
        else
        	error('Incorrect number of arguments.');
        end
        
        Random_Matrix(yr,:) = randi(nQ, 1, 52);
        % the bootstrapping value for each week of the synthetic record
        % is set at random based on the size of the 'perturbed
        % historical' record.  This is set before the real generating
        % loop because we want the same values for each streamflow
        % station (k loop below)
    end
    
    for k = 1:npoints
        for yr = 1:num_syn_years
            period = ceil(yr/period_length);
            % again, determine which perturbation period we are in
            if p == 0
                nQ = nQ_historical;
            elseif p > 0
                nQ = nQ_historical + ceil((p*nQ_historical))*n(period) + (nQ_historical - ceil((p*nQ_historical)))*m(period);
                % need to recalculate the size of the historical perturbed
                % record for the generation loop, same as above
            else
                error('Incorrect number of arguments.');
            end
    
            Q_matrix_int = Q_historical{k};%% getting the original streamflow matrix from streamflow station 'k'
            Q_matrix = Q_matrix_int;
            if p > 0
                temp = sort(Q_matrix_int);%%sort for biasing
                
                append = temp(1:ceil(p*nQ_historical),:); 
                % find lowest p% of values for each week
                
                append2 = temp((ceil(p*nQ_historical)+1):nQ_historical,:);
                % find the remaining streamflow values
                
                Q_matrix_int_int = vertcat(Q_matrix_int, repmat(append, n(period), 1));
                % append the lowest p% of streamflow values to the
                % end of the historical streamflow record, n times (value of n
                % changes depending on the perturbation period)
                
                Q_matrix = vertcat(Q_matrix_int_int,repmat(append2,m(period),1));
                % append the remaining highest (1-p)% to the historical +
                % lowflow record, m times (value of m changes depending on
                % the perturbation period)
            end            
            logQ = log(Q_matrix);   % log transform the perturbed record
            logQint = log(Q_matrix_int);    
            % log transform the original, unperturbed record 
            % (for means, SD, and autocorrelation calculations)

            Z = zeros(nQ, 52);

            for i=1:52
                weekly_mean(i) = mean(logQint(:,i));    
                % means of the original record
                weekly_stdev(i) = std(logQint(:,i));    
                % stdev of the original record
                Z(:,i) = (logQ(:,i) - weekly_mean(i)) / weekly_stdev(i);
                % whitened residual of the original record
            end
            for i=1:52
                Qs_uncorr(yr,i) = Z(Random_Matrix(yr,i), i);
                % bootstrap values for this year
            end
        end
        
        % this takes the residual record and shifts it by 26 weeks (so that
        % one square matrix starts each year one week 1, and another starts
        % each year at week 27
        Z_vector = reshape(Z',1,[]);    % makes it a 1D vector
        Z_shifted = reshape(Z_vector(27:(nQ*52-26)),52,[])';%shifts than reshapes to a 80X52 vector

        % The correlation matrices should use the historical Z's
        % (the "appended years" do not preserve correlation) 
        % correlation is calculated for the original and shifted values,
        % because complete autocorrelations can only be calculated
        % after a number of values have been observed in that year, so that
        % the original (starting at week 1) give us good autocorrelation in
        % weeks 27-52, and the shifted (starting at week 27) give us good
        % autocorrelation in weeks 1 - 26
        
        % Cholesky Decomposition - gives us the Upper Triangluar matricies
        U = chol(corr(Z(1:nQ_historical,:)));
        U_shifted = chol(corr(Z_shifted(1:nQ_historical-1,:)));
        
        Qs_uncorr_vector = reshape(Qs_uncorr(:,:)',1,[]);
        % uncorrelated, bootstrapped values as 1D vector
        Qs_uncorr_shifted(:,:) = reshape(Qs_uncorr_vector(27:(num_syn_years*52-26)),52,[])';
        % this shifts the uncorrelated, bootstrapped values so that both U
        % and U_shifted are matched with the same uncorrelated,
        % bootstrapped values (i.e. we can stitch the results together and
        % maintain a continuous correlation)
        
        % Takes the Upper Triangluar matricies and creates random matricies
        % that take on the autocorrelation structure of the original Z and
        % Z_shifted matricies
        Qs_corr(:,:) = Qs_uncorr(:,:)*U;
        Qs_corr_shifted(:,:) = Qs_uncorr_shifted(:,:)*U_shifted;

        % Stitch together the matricies so that we are only using the
        % values with a complete picture of autocorrelation (the back half
        % of each matrix)
        Qs_log(:,1:26) = Qs_corr_shifted(:,27:52);
        Qs_log(:,27:52) = Qs_corr(2:num_syn_years, 27:52);
        
        % add back in the mean, std, and log transformation
        for year = 1:(num_syn_years-1)
            for i=1:52
                Qsk(year,i) = exp(Qs_log(year,i)*weekly_stdev(i) + weekly_mean(i));
            end
        end
        % do for each streamflow station
        Qs{k} = Qsk;
    end
end
