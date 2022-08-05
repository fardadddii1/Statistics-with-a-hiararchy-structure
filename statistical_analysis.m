clear all
close all

n_nodes = [25 50 75 100 150]; %Number of cross sections

n_samples = 35; %number of samples on each cross section
n_regions = 5; % number of regions

% cross_info1 = [cross_info, shape_variance', perimeter, circle_perimeter', cross_area',circle_area', maxjor_ax_length', minor_ax_length', eccentricity', compactness'];

% Features3 = array2table(cross_info1,...
%     'VariableNames',{'dmax', 'dmin', 'dmean', 'dstd', 'shape variance', 'perimeter', 'template circle perimeter', 'cross section area','template circle area', 'maxjor ax length', 'minor ax length', 'eccentricity (minor/major)', 'compactness'});

save_files = 'C:\Users\smousavi\Desktop\fardad\28-vessel-morphometric\results\final\for_statistics';

output_path = strcat(save_files,'\feat_results');
if ~exist(output_path);
    mkdir(output_path);
end

cross_section = 'cross_section_';

for ui = 1: length(n_nodes)
    %     cross_section_names = [];
    nn = n_nodes(ui);
    
    node_codes=cell(nn,13);
    
    Failed_cases_ad = append(save_files, '\',string(nn),'_failed_featutres_stacked_by_cases_first.csv');
    
    Success_cases_ad = append(save_files, '\',string(nn),'_successfull_featutres_stacked_by_cases_first.csv');
    
    Failed_cases = readtable(Failed_cases_ad);
    
    Success_cases = readtable(Success_cases_ad);
    
    cross_section_names = append(cross_section,string(linspace(1,nn,nn)));
    
    for kl=1:nn
        
        succ_idx = Success_cases{:,'section_number'}==cross_section_names(kl);
        fail_idx = Failed_cases{:,'section_number'}==cross_section_names(kl);
        
        succ_feats(:,:, kl) = table2array(Success_cases(succ_idx,3:end));
        fail_feats(:,:, kl) = table2array(Failed_cases(fail_idx,3:end));
        
        feat_stacked = [succ_feats(:,:, kl); fail_feats(:,:, kl)];
        
        group_stacked = [ones(size(succ_feats,1),1);zeros(size(fail_feats,1),1)];
        group_stacked_1_2 = [ones(size(succ_feats,1),1);zeros(size(fail_feats,1),1)+2];
        for ty=1:size(succ_feats,2)
            
            node_codes(kl,ty) = {strcat(num2str(kl),'_',num2str(ty))};
            
            %             Normality_test_sucess(:,:,ty) = normalitytest (succ_feats(:,ty,kl)');
            %             Normality_test_fail(:,:,ty) = normalitytest (fail_feats(:,ty,kl)');
            [h_suc(ty),p_suc(ty)] = adtest(succ_feats(:,ty, kl), 'Alpha', 0.05);
            [h_fail(ty),p_fai(ty)] = adtest(fail_feats(:,ty, kl), 'Alpha', 0.05);
            
            if h_suc(ty) == 1 || h_fail(ty) ==1
                
                [p_equ(ty), stats_equ{ty}] = vartestn(feat_stacked(:,ty), group_stacked, 'display','off', 'TestType', 'LeveneQuadratic');
                
                if p_equ(ty)<0.05
                    
                    [p_sig(ty), h_sig(ty),stats_sig{ty}] = ranksum(succ_feats(:,ty, kl),fail_feats(:,ty, kl));
                    
                else
                    
                    [p_sig(ty), h_sig(ty),stats_sig{ty}] = ranksum(succ_feats(:,ty, kl),fail_feats(:,ty, kl));
                    
                end
                
            else
                
                [p_equ(ty), stats_equ{ty}] = vartestn(feat_stacked(:,ty), group_stacked, 'display','off', 'TestType', 'Bartlett');
                
                if p_equ(ty)<0.05
                    
                    p_sig(ty) = Wtest([feat_stacked(:,ty), group_stacked_1_2]);
                    
                else
                    
                    [h_sig(ty),p_sig(ty),chi, stats_sig{ty}] = ttest2(succ_feats(:,ty, kl),fail_feats(:,ty, kl));
                    
                end
                
            end
            
            
        end
        
        p_final(kl,:) = p_sig;
        
    end
    
    results = [node_codes(p_final<0.05),num2cell(p_final(p_final<0.05))];
    results = array2table(results);
    feat_names = results{:,1};
    writetable(results,strcat(output_path,'\Sig_Feats_',num2str(nn),'_nodes.csv'));
    
    succ_feat_subset = [];
    fail_feat_subset = [];
    
    succ_temp = Success_cases{:,3:end};
    fail_temp = Failed_cases{:,3:end};
    succ_labels = [];
    fail_labels = [];
    for j = 1:length(results{:,1});
        feat_id = strsplit(results{j,1}{1},'_');
        cs_idx = string(strcat(cross_section,feat_id{1}));
        feat_idx = str2num(feat_id{2});
        
        succ_rows = Success_cases{:,'section_number'}==cs_idx;
        fail_rows = Failed_cases{:,'section_number'}==cs_idx;
        succ_cols = succ_temp(succ_rows,feat_idx);
        fail_cols = fail_temp(fail_rows,feat_idx);
        succ_labels = Success_cases{Success_cases{:,'section_number'}==cs_idx,'case_ID'};
        fail_labels = Failed_cases{Failed_cases{:,'section_number'}==cs_idx,'case_ID'};
        succ_feat_subset = [succ_feat_subset,succ_cols];
        fail_feat_subset = [fail_feat_subset,fail_cols];
    end
    
    results = array2table(succ_feat_subset,'RowNames',succ_labels,'VariableNames',feat_names);
    writetable(results,strcat(output_path,'\Sig_Feat_Vals_Success_',num2str(nn),'_nodes.csv'),'WriteRowNames',true,...,
        'WriteVariableNames',true);
    
    results = array2table(fail_feat_subset,'RowNames',fail_labels,'VariableNames',feat_names);
    writetable(results,strcat(output_path,'\Sig_Feat_Vals_Failure_',num2str(nn),'_nodes.csv'),'WriteRowNames',true,...,
        'WriteVariableNames',true);
    
end

