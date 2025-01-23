%% Read data
% used labels:
% HIV pos y==1
% HIV neg y==0

close all
clear
filename = 'Data/BioPlex_PCR DEID Data - Complete v3.xlsx';
filename2 = 'Data/matched_test_results_20160106-20171031.csv';
filename3 = 'Data/matched_test_results_20171101-20181231.csv';
[~,~,all_1] = xlsread(filename,'BioPlex Final Disp');
[~,~,tmp1] = xlsread(filename2);
[~,~,tmp2] = xlsread(filename3);
all_2 = [tmp1(2:end, :) ;tmp2(2:end, :)];
clear filename filename2 filename3 tmp1 tmp2

%% Process data
X1 = all_1(2:end ,4:7);
Y1 = all_1(2:end,13);
ID = all_1(2:end ,1);
assay_results = all_1(2:end,14);
rna_stat = all_1(2:end,11);
ID_1_org = ID;
% remove entries with 'IND' labels in final disposition
rows = cellfun(@(x) strcmp('IND',x),Y1);
X1(rows,:) = [];
ID(rows) = [];
assay_results(rows) = [];
% 1 positive, 0 negative
Y1 = double(cellfun(@(x) strcmp('+',x),Y1)); % negative == 0 positive == 1
Y1(rows) = []; %rows where there were IND still in 'rows', DO NOT CHANGE ORDER
assay_results = double(cellfun(@(x) strcmp('+',x),assay_results)); % negative ==0 positive ==1

% remove chars (indeterminate values)
tmp = cellfun(@(x) ischar(x),X1);
[rows,~] = find(tmp);
X1(rows,:) = [];
Y1(rows) = [];
ID(rows) = [];
assay_results(rows) = [];
clear tmp rows 

X1 = cell2mat(X1);

% %%%% Outliers HIV2 Ab
% id = X1(:,4)>200;
% Y1(id) = nan;
% X1(id,:) = nan;
% ID(id) = [];
% rna_stat(id) = nan;
% essay_result(id) = nan;


% combining data
combine = true;
perc = 0.25 ; % percentage of data to take randomily from data other than
% Bioplex file. We remove entries with values > 1
if combine
    ID2 = all_2(2:end ,1);
    X2  = all_2(2:end, 3:6);
    % remove IDs existing in Bioplex files from the two other files
    rows = ismember(ID2,ID_1_org);
    X2(rows,:) = [];
    ID2(rows,:) = [];
    %     % rermove chars (indeterminate values), None were found, commenting it
    %     tmp = cellfun(@(x) ischar(x),X2); %very slow
    %     [rows,~] = find(tmp);
    %     X2(rows,:) = [];
    %     ID2(rows,:) = [];
    X2 = cell2mat(X2);
    X2_orig = X2;
    
    IDX = randi(length(X2), ceil(perc*length(X2)),1);
    ID2 = ID2(IDX);
    X2 = X2(IDX,:);
    % above 1 on HIV_Ag_Ab considered positive??
    %Y2 = X2(:,1) >= 1;
    % removing values bigger than 1 the threshold
    X2(X2(:,1)>1,:) = [];
    Y2 = zeros(length(X2),1);
    second_neg_id = length(X1)+1:length(X1)+length(X2);
    
    X = [X1; X2];
    Y = [Y1; Y2];
    %ID = [ID;ID2];
    
end
pos_id = Y==1;
neg_id = Y==0;
clear ID ID2 IDX id combine perc rows tmp thresh
Features = [{'HIV Ag Ab'},{'HIV 1 Ab'},{'HIV 1 Ag'},{'HIV 2 Ab'}];

%% classify with 3 params on HIV_1_Ab, HIV_1_Ag, HIV_2_Ab
rng default
features = "HIV 1 Ab, HIV 1 Ag, HIV 2 Ab";
clear X_classify Y_classify
% training_labels_type: changes the used output labels in the classification
% "ground truth = 1" or "assay results = 2" or both = 3,
% If both is chosen as Y, predicted labels and P_vals etc. will have 2 columns
% making predicted labels and entropy vals first column: essay result
% labels, second column: gound truth labels
training_labels_type = 2;
numOfFolds = 10;

X_classify = X;
if training_labels_type==1
    Y_classify = Y;
    n = 1;
elseif training_labels_type==2
    Y_classify = [assay_results; zeros(length(X2),1)];
    n = 1;
else
    n = 2;
    Y_classify(:,1) = [assay_results; zeros(length(X2),1)];
    Y_classify(:,2) = Y;
end

% --- Assay classifier vs final disposition confusion table

conf_sp = confusionmat(Y(1:length(X1)), assay_results);

predicted_labels = zeros( length(X_classify), n);
scores = zeros( length(X_classify), n);
accuracy = zeros(1,n); sensitivity = zeros(1,n); specificity = zeros(1,n);
false_positive_id=zeros(length(Y_classify),n);false_negative_id=zeros(length(Y_classify),n);
total_conf = zeros(2,2,2);
special_conf = zeros(2,2);%assay classifier vs final disposition

% Classification and Cross-validation 
for i=1:n
    c = cvpartition(Y_classify(:,i),'k',numOfFolds);
    for j=1:numOfFolds
        train_data = X_classify(c.training(j),2:4);
        test_data = X_classify(c.test(j),2:4);
        train_labels = Y_classify(c.training(j),i);
        test_labels = Y_classify(c.test(j),i);
        %classifier
        svm = fitcsvm(train_data,train_labels,'KernelFunction','rbf',...
    'Standardize',false);
        [predicted_labels_tmp,tmp_scores] = predict(svm, test_data);
%         p = 1./(1+exp(-abs(tmp_scores(:,1))));
%         p_vals(c.test(j),i) = p;
        scores(c.test(j),i) = tmp_scores(:,1);
        predicted_labels(c.test(j),i)= predicted_labels_tmp;
        
        conf = confusionmat(test_labels, predicted_labels_tmp);
        total_conf(:,:,i) = total_conf(:,:,i) + conf;     
        
    end
    
    %mean accuracy sensitivity and specificity
    accuracy(i) = (total_conf(1,1,i)+total_conf(2,2,i))*100/sum(total_conf(:,:,i), 'All');
    sums = sum(total_conf(:,:,i), 2);
    specificity(i) = total_conf(2,2,i)*100/sums(2);
    sensitivity(i) = total_conf(1,1,i)*100/sums(1); % Accuracy of detection for negative cases

    %Save "prediction wise" false positive and false negative
    
    false_positive_id(:,i) = predicted_labels(:,i) ~= Y_classify(:,i) & predicted_labels(:,i)==1;
    false_negative_id(:,i) = predicted_labels(:,i) ~= Y_classify(:,i) & predicted_labels(:,i)==0;
end

p_vals = 1./(1+exp(-(abs(scores)-mean(scores,'omitnan'))*10));
entropy_ = @(x) -1*(x*log2(x)+(1-x)*log2(1-x));
entropy_vals = arrayfun(entropy_, p_vals);

% ---1)--- assay classifier vs final disposition confusion table (first file)
% ---2)--- assay classifier vs final disposition confusion table (first file)

% first_file_conf = zeros(2,2,2); % 1 assay classifier 2 clilnical classifier vs final disposition
% first_file_conf(:,:,1) = confusionmat(Y(1:length(X1)),predicted_labels(1:length(X1),1));
% first_file_conf(:,:,2) = confusionmat(Y(1:length(X1)),predicted_labels(1:length(X1),2));
% --- assay results vs final disposition confusion table
assay_results_conf = confusionmat(Y(1:length(X1)), assay_results);
% --- assay classifier vs final disposition confusion table (first file)
% scatter3(X_(:,1),X_(:,2),p_vals(:,1))
clear X_classify Y_classify p idx predicted_labels_tmp conf z train_labels test_labels train_data test_data  j sums c

%% --- USING ONLY positive assay entries.. 
% For first file only, with Assay definition of false cases

close all
entropy_ = @(x) -1*(x*log2(x)+(1-x)*log2(1-x));
entropy_vals = arrayfun(entropy_, p_vals);
% --entropy vals used to be really for entropies, butnow it is a proxy for
% features used for classification--
% entropy_vals_tmp = entropy_vals(1:sum(idx),1); % these come from first file AND do not have NaN entries
scores_vals_tmp = scores(1:length(X1),1);

true_pos = Y1==1 & assay_results == Y1;
false_pos = Y1==0 & assay_results ~= Y1;
true_neg = Y1==0 & assay_results == Y1;
false_neg = Y1==1 & assay_results ~= Y1;

h=figure;
set(h,'color','w')
tiledlayout(2,1);

% Assay Positive cases, Assay classifier entropies
ax1 = nexttile;
histogram(scores_vals_tmp(true_pos),20, 'Normalization', 'Probability','FaceColor','#009F28', 'FaceAlpha', 0.5);
mean_ = num2str(mean(scores_vals_tmp(true_pos), 'omitnan'));
std_ = num2str(std(scores_vals_tmp(true_pos), 'omitnan'));
txt = {"True positive:",strcat("mean=", mean_,", std=", std_)};
text(0.1,0.8,txt,'Units','Normalized','Color','#009F28')
hold on
scatter(scores_vals_tmp(true_pos), 0.4+0.35*rand(1,length(scores_vals_tmp(true_pos))), 'MarkerEdgeColor','#009F28', 'MarkerFaceColor','#009F28')

h=histogram(scores_vals_tmp(false_pos ),30, 'Normalization', 'Probability','FaceColor','r', 'FaceAlpha', 0.5);
mean_ = num2str(mean(scores_vals_tmp(false_pos),'omitnan'));
std_ = num2str(std(scores_vals_tmp(false_pos),'omitnan'));
txt = {"False positive:",strcat("mean=", mean_,", std=", std_)};
text(0.1,0.6,txt,'Units','Normalized','Color','r')

% predicted labels of Assay classifier
predicted_labels_tmp = predicted_labels(1:length(X1),1);
% false positive cases from assay classifier that are Falsely
% predicted by model Assay classifier
idx_ = false_pos & (predicted_labels_tmp ~= Y1);
scatter(scores_vals_tmp(idx_), 0.2+0.2*rand(1,length(scores_vals_tmp(idx_))), 'MarkerEdgeColor','r', 'MarkerFaceColor','r','Marker','o')
% false positive cases from training on essay results, and also CORRECTLY
% predicted by model from training on assay results
idx_ = false_pos & (predicted_labels_tmp == Y1);
scatter(scores_vals_tmp(idx_), 0.2+0.35*rand(1,length(scores_vals_tmp(idx_))), 'MarkerEdgeColor','r', 'MarkerFaceColor','r','Marker','^','MarkerFaceAlpha',0.2)
legend([{'Ture + distribution '},{'True +'},{'False + distribution'},{'False + falsely predicted by Assay Classifier'},{'False + correctly predicted by Assay Classifier'}])
title({'Positive Assay Results cases scores distributions from Assay Classifier', 'First file'})
xlim([-inf,inf])
xlabel("Score", 'FontSize', 13)
ylabel("Probability", 'FontSize', 13)
set(gca,'FontSize',13)


% Assay Positive cases, Clinical classifier entropies
% entropy_vals_tmp = entropy_vals(1:sum(idx),2); % these come from first file AND do not have NaN entries
% entropy_vals_tmp = scores(1:sum(idx),2);
scores_vals_tmp = scores(1:length(X1),2);
ax2 = nexttile;
histogram(scores_vals_tmp(true_pos),20, 'Normalization', 'Probability','FaceColor','#009F28', 'FaceAlpha', 0.5);
mean_ = num2str(mean(scores_vals_tmp(true_pos), 'omitnan'));
std_ = num2str(std(scores_vals_tmp(true_pos), 'omitnan'));
txt = {"True positive:",strcat("mean=", mean_,", std=", std_)};
text(0.1,0.4,txt,'Units','Normalized','Color','#009F28')
hold on
scatter(scores_vals_tmp(true_pos), 0.4+0.35*rand(1,length(scores_vals_tmp(true_pos))), 'MarkerEdgeColor','#009F28', 'MarkerFaceColor','#009F28')
hold on
h_=histogram(scores_vals_tmp(false_pos),20, 'Normalization', 'Probability','FaceColor','r', 'FaceAlpha', 0.5);
mean_ = num2str(mean(scores_vals_tmp(false_pos),'omitnan'));
std_ = num2str(std(scores_vals_tmp(false_pos),'omitnan'));
txt = {"False positive:",strcat("mean=", mean_,", std=", std_)};
text(0.1,0.2,txt,'Units','Normalized','Color','r')
hold on
% predicted labels from training on final disposition
predicted_labels_tmp = predicted_labels(1:length(X1),2);
% false positive cases from training on essay results, and also Falsely
% predicted by model from training on final disposition
idx_ = false_pos & (predicted_labels_tmp ~= Y1);
scatter(scores_vals_tmp(idx_), 0.2+0.2*rand(1,length(scores_vals_tmp(idx_))), 'MarkerEdgeColor','r', 'MarkerFaceColor','r','Marker','o')
hold on
% false positive cases from training on essay results, and also CORRECTLY
% predicted by model from training on final disposition
idx_ = false_pos & (predicted_labels_tmp == Y1);
scatter(scores_vals_tmp(idx_), 0.45+0.35*rand(1,length(scores_vals_tmp(idx_))), 'MarkerEdgeColor','r', 'MarkerFaceColor','r','Marker','^','MarkerFaceAlpha',0.2)


legend([{'Ture + distribution '},{'True +'},{'False + distribution'},{'False + falsely predicted by Clinical Classifier'},{'False + correctly predicted by Clinical Classifier'}])
title({'Positive Assay Results cases scores distributions from Clinical Classifier', 'First file'})
xlim([-inf,inf])
xlabel("Score")
ylabel("Probability")
set(gca,'FontSize',13)

%% Using only Assay positive results
% predict if the positive cases of the essay results are false positive
% or true positive
close all

% taking cases predicted positive by assay results 
% (Assay results: where if max(features)>1 it is considered a positivive
% case, otherwise negative)

assay_positive_idx = assay_results == 1;
Y_tmp = Y1(assay_positive_idx)==0;

PC_ASSAY_CHOICE = 1; %1 2 or 3 for PC1 PC2 Assay meas. respectively
if any(PC_ASSAY_CHOICE == [1, 2])
    X_tmp_ = scores(1:length(Y1),:);
    % Doing PCA for assay positive cases
    X_tmp_ = X_tmp_(assay_positive_idx,:);
    [~,Xm_tp_,~] = pca(X_tmp_);
    % HERE choose All ":", PC1 "1", PC2 "2"
    X_tmp = X_tmp_(:,PC_ASSAY_CHOICE);
    if PC_ASSAY_CHOICE==1
        Penalty = 1.5;
    else
        Penalty = 4;
    end
elseif PC_ASSAY_CHOICE == 3
    X_tmp = X1(:,2:4);
    X_tmp = X_tmp(assay_positive_idx,:);
    Penalty = 1;
else
    X_tmp_=0;%error
end

total_conf_2 = zeros(2,2); % total confusion matrix : sum of confusion
% matrices from each fold
numOfFolds = 10;
c = cvpartition(Y_tmp,'k',numOfFolds);
cost_vals = linspace(1,4,50)';
accuracy_2 = zeros(length(cost_vals),1); sensitivity_2 = zeros(length(cost_vals),1); specificity_2 = zeros(length(cost_vals),1);
%%
for ii= 1:length(cost_vals)
    total_conf_2 = zeros(2,2);
    for j=1:numOfFolds
        train_data = X_tmp(c.training(j));
        test_data = X_tmp(c.test(j));
        train_labels = Y_tmp(c.training(j));
        test_labels = Y_tmp(c.test(j));
        %     svm = fitcsvm(train_data,train_labels,'KernelFunction','polynomial', 'Standardize',true,'Cost',[0 1; 3.1 0],'PolynomialOrder',5);%, 'Standardize',true);
        svm = fitcsvm(train_data,train_labels,'KernelFunction','rbf', 'Standardize',true,'Cost',[0 1; cost_vals(ii) 0]);%, 'BoxConstraint',4.3507, 'KernelScale',0.67278);
        predicted_labels_entropy_classifier = predict(svm, test_data);
        [conf,ord] = confusionmat(test_labels, predicted_labels_entropy_classifier);
        total_conf_2 = total_conf_2 + conf;
    end
    
    %mean accuracy sensitivity and specificity
    accuracy_2(ii) = (total_conf_2(1,1)+total_conf_2(2,2))*100/sum(total_conf_2, 'All');
    sums = sum(total_conf_2, 2);
    sensitivity_2(ii) = total_conf_2(2,2)*100/sums(2); %accuracy of False positive prediction
    specificity_2(ii) = total_conf_2(1,1)*100/sums(1);
end

figure('color','w')
plot(cost_vals, [accuracy_2, sensitivity_2, specificity_2])
xlabel("Sensitivity penalty factor (to improve sensitivity)")
legend(["Accuracy","Sensitivity","Specificity"])
set(gca,'FontSize',13)
% -scatter
figure('color','w')
scatter(X_tmp_(~Y_tmp,1),X_tmp_(~Y_tmp,2), 'MarkerFaceColor','g' )
hold on
scatter(X_tmp_(logical(Y_tmp),1),X_tmp_(logical(Y_tmp),2), 'MarkerFaceColor','r')
legend([{'True +'},{'False +'}])
% xlabel('entropy from essay results classifier')
% ylabel('entropy from final disposition classifier')
xlabel('PC1')
ylabel('PC2')
set(gca,'FontSize',13)
%draw pc1
figure('color','w')
scatter(X_tmp_(~Y_tmp,1), zeros(length(X_tmp_(~Y_tmp,1)),1), 'MarkerFaceColor','g' )
hold on
scatter(X_tmp_(logical(Y_tmp),1),zeros(length(X_tmp_(logical(Y_tmp),1)),1), 'MarkerFaceColor','r')
legend([{'True +'},{'False +'}])
xlabel('PC1')
set(gca,'FontSize',13)
%draw pc2
figure('color','w')
scatter(X_tmp_(~Y_tmp,2), zeros(length(X_tmp_(~Y_tmp,2)),1), 'MarkerFaceColor','g' )
hold on
scatter(X_tmp_(logical(Y_tmp),2),zeros(length(X_tmp_(logical(Y_tmp),2)),1), 'MarkerFaceColor','r')
legend([{'True +'},{'False +'}])
xlabel('PC2')
set(gca,'FontSize',13)
%%
% rng default
rep_count = 10;
accuracy_2 = zeros(rep_count,1); sensitivity_2 = zeros(rep_count,1); specificity_2 = zeros(rep_count,1);
opts = struct('Optimizer','bayesopt','ShowPlots',false,...
    'AcquisitionFunctionName','expected-improvement-plus','Verbose',0);
accuracy_fold = struct();


for ii=1:rep_count
        rng shuffle
        disp({'Count ',ii})
        total_conf_2 = zeros(2,2);
        c = cvpartition(Y_tmp,'k',numOfFolds);
    for j=1:numOfFolds
        disp({'Fold',j})
        train_data = X_tmp(c.training(j));
        test_data = X_tmp(c.test(j));
        train_labels = Y_tmp(c.training(j));
        test_labels = Y_tmp(c.test(j));
    %     svm = fitcsvm(train_data,train_labels,'KernelFunction','polynomial', 'Standardize',true,'Cost',[0 1; 3.1 0],'PolynomialOrder',5);%, 'Standardize',true);
        svm = fitcsvm(train_data,train_labels,'KernelFunction','rbf', 'Standardize',true,'Cost',[0 1; Penalty 0],'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts);
        predicted_labels_ = predict(svm, test_data);
        [conf,ord] = confusionmat(test_labels, predicted_labels_);
        total_conf_2 = total_conf_2 + conf; 
        accuracy_fold.(strcat(['Rep_',num2str(ii)])).Accuracy(j) = (conf(1,1)+conf(2,2))*100/sum(conf, 'All');
        sums = sum(conf, 2);
        accuracy_fold.(strcat(['Rep_',num2str(ii)])).Sensitivity(j) = conf(2,2)*100/sums(2); %accuracy of False positive prediction
        accuracy_fold.(strcat(['Rep_',num2str(ii)])).Specificity(j) = conf(1,1)*100/sums(1);

    end

    %mean accuracy sensitivity and specificity
    accuracy_2(ii) = (total_conf_2(1,1)+total_conf_2(2,2))*100/sum(total_conf_2, 'All');
    sums = sum(total_conf_2, 2);
    sensitivity_2(ii) = total_conf_2(2,2)*100/sums(2); %accuracy of False positive prediction
    specificity_2(ii) = total_conf_2(1,1)*100/sums(1);
end

means = mean([accuracy_2, sensitivity_2, specificity_2]);
stds = std([accuracy_2, sensitivity_2, specificity_2]);

