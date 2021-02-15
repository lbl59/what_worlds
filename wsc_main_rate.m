clc; clear all;

% all input inflows are in million gallons/week
inflow_dir = 'historical-data';
inflow_files = {'trainingLittleRiverRaleighInflow';'trainingOWASAInflow';'trainingClaytonGageInflow';'trainingCrabtreeCreekInflow';'trainingFallsLakeInflow';
    'trainingJordanLakeInflow';'trainingLakeWBInflow';'trainingLillingtonInflow';'trainingLittleRiverInflow';'trainingMichieInflow'};               

num_syn_realizations = 1000;    % number of synthetic realizations
num_syn_years = 1;     % number of synthetic years
for k=1:10
    
	Qh{k} = load([inflow_dir '/' inflow_files{k} '.csv']);    
	output_dyn{k} = zeros(num_syn_realizations, num_syn_years*52);
    output_stat{k} = zeros(num_syn_realizations, num_syn_years*52);
    
end

% generate realizations
p_dyn = 0.25;    % low X% of streamflows
p_stat = 0;
n = [1 1 1 1 1 1];     % how much more likely are those low flows? (1 = same as historical)
m = [1 1 8 12 20 20];     % how much more likely are the remaining flows? (i.e. wettest (1 - X)%)

for r=1:num_syn_realizations
	disp(r);
	Qs_dyn = stress_dynamic(Qh, num_syn_years, p_dyn, n, m);
    Qs_stat = stress_dynamic(Qh, num_syn_years, p_stat, n, m);
    for k=1:10
        output_dyn{k}(r,:) = reshape(Qs_dyn{k}',1,[]);
        output_stat{k}(r,:) = reshape(Qs_stat{k}',1,[]);
	end

end
        
for k=1:10
	dlmwrite(['synthetic-data-dyn/' inflow_files{k} '_SYN01.csv'], output_dyn{k});
    dlmwrite(['synthetic-data-stat/' inflow_files{k} '_SYN01.csv'], output_stat{k});
end
        
    