%% Classification using a threshold
% Here we test the performance of the classification only by using a
% threshold, and show ROC curve
close all
clear
filename = 'Data/BioPlex_PCR DEID Data - Complete v3.xlsx';
filename2 = 'Data/matched_test_results_20160106-20171031.csv';
filename3 = 'Data/matched_test_results_20171101-20181231.csv';
[~,~,all_1] = xlsread(filename,'BioPlex Final Disp');
[~,~,tmp1] = xlsread(filename2);
[~,~,tmp2] = xlsread(filename3);
varnames = tmp1(1, :);
all_2 = [tmp1(2:end, :) ;tmp2(2:end, :)];
clear filename filename2 filename3

%% prepare data
X1 = all_1(2:end ,4); %only assay
Y1 = all_1(2:end,13); % HIV Final disposition
ID = all_1(2:end ,1);
rna_stat = all_1(2:end,11);
% remove entries with 'IND' labels in final disposition
rows1 = cellfun(@(x) strcmp('IND',x),Y1);
X1(rows1,:) = [];
Y1(rows1) = []; %rows where there were IND still in 'rows', DO NOT CHANGE ORDER
ID_1_org = ID;
ID(rows1) = [];
% 1 positive, 0 negative
Y1 = double(cellfun(@(x) strcmp('+',x),Y1)); % negative ==0 positive ==1

% % remove chars (indeterminate values) !! Don't
% tmp = cellfun(@(x) ischar(x),X1);
% [rows,~] = find(tmp);
% X1(rows,:) = [];
% Y1(rows) = [];

X1 = cell2mat(X1);
% combining data
combine = true;
if combine
    all_2_mod = all_2;
    ID2 = all_2(: ,1);
    X2  = all_2(:, 3); %only assay
    
    % remove IDs existing in Bioplex files from the two other files
    rows = ismember(ID2,ID_1_org);
    X2(rows,:) = [];
    ID2(rows,:) = [];
    all_2_mod(rows,:) = [];
    %     % rermove chars (indeterminate values),NOTE None were found, commenting it
    %     tmp = cellfun(@(x) ischar(x),X2); %very slow
    %     [rows,~] = find(tmp);
    %     X2(rows,:) = [];
    %     ID2(rows,:) = [];
    X2 = cell2mat(X2);
    % above 1 on HIV_Ag_Ab considered positive??
    %Y2 = X2(:,1) >= 1;
    
%     odd_cases = all_2_mod(X2(:,1) >= 1,:);
%     odd_ids = odd_cases(:,1);

    % removing values bigger than 1 the threshold... SHOULD I?? 12/14/2021
%     X2(X2(:,1)>1,:) = [];

    % remove specific IDs
    rows = ismember(ID2,["xxxx", "xxxx", "xxxx"]); % Removed IDs for privacy (online)
    X2(rows,:) = [];
    ID2(rows,:) = [];
    all_2_mod(rows,:) = [];
   
    Y2 = zeros(length(X2),1);
    
    X = [X1; X2];
    Y = [Y1; Y2];
    %ID = [ID;ID2];
else
    X = X1;
    Y = Y1;
end
clear tmp rows non_nan_rna_id 
% T = cell2table(odd_cases, 'VariableNames', varnames);
% writetable(T, '.\Threshold test and ROC\odd_cases_2_3rd_files.csv')
%%
assay = X;
Threshold = linspace(0,ceil(max(assay))+1,ceil(max(assay))+100)';

Accuracy = zeros(size(Threshold));
Sensitivity = zeros(size(Threshold));
Specificity = zeros(size(Threshold));

for k = 1:length(Threshold)
    assay_res = zeros(size(assay));
    assay_res(assay>=Threshold(k)) = 1;
    
    [conf,~] = confusionmat(Y, assay_res);
    Accuracy(k) = (conf(1,1)+conf(2,2))*100/sum(conf, 'All');
    sums = sum(conf, 2);
    Sensitivity(k) = conf(2,2)*100/sums(2); %accuracy of False positive prediction
    Specificity(k) = conf(1,1)*100/sums(1);
end
% a=table(Threshold,Accuracy,Sensitivity,Specificity);
% writetable(a,'assay_threshold_test.xls')
%%
% close all
figure('color','w')
plot(Threshold,[Accuracy,Sensitivity,Specificity],'linewidth',1.1)
legend(["Accuracy","Sensitivity","Specificity"],'location','southeast')
xlabel("Threshold")
ylabel("Percentage %")
set(gca,'FontSize',14)
xlim([0,100])
ylim([90,100])
%%
 close all
J_stat = Sensitivity + Specificity - 100;
FPR = 100-Specificity;
[ ~, idx] =  max(J_stat);
best_threshold = Threshold(idx);
figure('color','w')
plot(flip(FPR), flip(Sensitivity),'linewidth',1.1)
hold on 
plot(FPR(idx), Sensitivity(idx), 'x', 'linewidth',2,'MarkerSize',10)
text(FPR(idx)+0.1, Sensitivity(idx)-0.2, join(['Threshold = ', num2str(best_threshold)]))
% legend(["Accuracy","Sensitivity","Specificity"],'location','southeast')
xlabel("FPR %")
ylabel("TPR %")
set(gca,'FontSize',14)
xlim([0,4])
ylim([95,100])