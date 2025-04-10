%%   classify  - 10 fold CV
clear;clc
result_path='/media/useradmin/databoxone/YCX/20230919_MovieGraphtheory/result/FCA';
mkdir(result_path)
cd(result_path)
%%%%% construct data
sup_path='/media/useradmin/databoxone/YCX/20230919_MovieGraphtheory/result//FIX-Denoised_200';
run_dir=dir(sup_path);
run_dir(1:2)=[];
label=1;
FC_all=zeros(169,15,10,200,200);
clips_num=[4,3,4,3];
tic
for irun=5:8
    clip_dir=dir(fullfile(sup_path,run_dir(irun).name,'clips_originvalue'));
    clip_dir(1:2)=[];
    for iclip=1:clips_num(irun-4)
        disp(strcat('loading...RunNum_',num2str(irun),'...Clip_',num2str(iclip)));
        sub_dir=dir(fullfile(clip_dir(iclip).folder,clip_dir(iclip).name,'GraphTheory_sparsity'));
        sub_dir(1:2)=[];
        for isub=1:length(sub_dir)
            load(fullfile(sub_dir(isub).folder,sub_dir(isub).name));
            for ispar=1:length(fc_sparsity)
                FC_all(isub,label,ispar,:,:)=cell2mat(fc_sparsity(ispar));
            end
        end
        label=label+1;
        toc
    end
end
% validation clip
for irun=5:8
    clip_dir=dir(fullfile(sup_path,run_dir(irun).name,'clips_originvalue'));
    clip_dir(1:2)=[];
    for iclip=clips_num(irun-4)+1
        disp(strcat('loading...RunNum_',num2str(irun),'...Clip_',num2str(iclip)));
        sub_dir=dir(fullfile(clip_dir(iclip).folder,clip_dir(iclip).name,'GraphTheory_sparsity'));
        sub_dir(1:2)=[];
        for isub=1:length(sub_dir)
            load(fullfile(sub_dir(isub).folder,sub_dir(isub).name));
            for ispar=1:length(fc_sparsity)
                vc(isub,ispar,irun-4,:,:)=cell2mat(fc_sparsity(ispar));
            end
        end
        toc
    end
end
FC_all(:,label,:,:,:)=squeeze(mean(vc,3));
%
%%%%%%%  Classify
N=200;   % number of rois
nmEdges=(N*(N-1))/2;   % number of edges
nmSbj=169;   % number of subjects
nmScans=15;   % number of scans
n_repeats = 100;
folds = 10;
%%% SVM
all_accuracies_acsc=zeros(3,n_repeats,folds,7);      % Preallocate accuracy and F1
all_f1_scores_acsc=zeros(3,n_repeats,folds,7);
all_accuracies_full=zeros(3,n_repeats,folds);
all_f1_scores_full=zeros(3,n_repeats,folds);
for ispar=4:6
    % Perform 10 fold cross-validation
    % Repeat 100 times
    scores_ACSC_all=zeros(nmEdges,n_repeats,folds);
    parfor r = 1:n_repeats
        % Create a new partition for each repetition
        cv = cvpartition(nmSbj, 'KFold', folds);
        
        acc_fold_acsc = zeros(folds,7);
        f1_fold_acsc = zeros(folds,7);
        acc_fold_full = zeros(folds,1);
        f1_fold_full = zeros(folds,1);
        for i = 1:folds
            timeNow=clock;
                fprintf('-SparNum  %d -- Repetition %d/%d --- fold %d/%d \n', ispar,r, n_repeats,i,folds);
                disp(strcat(num2str(timeNow(1)),'-',num2str(timeNow(2)),'-',num2str(timeNow(3)),';;;',...
                    num2str(timeNow(4)),':',num2str(timeNow(5)),':',num2str(timeNow(6))));

            trainIdx = cv.training(i);
            testIdx = cv.test(i);
            Ytrain=categorical(reshape(repmat([1:15],sum(trainIdx),1),[],1));
            Ytest=categorical(reshape(repmat([1:15],sum(testIdx),1),[],1));
            
            % Split data into training and test sets
            FCtrain=squeeze(FC_all(trainIdx,:,ispar,:,:));
            FCtest=squeeze(FC_all(testIdx,:,ispar,:,:));
            
            % ACSC feature selection
            %%% choose features by ACSC
            flatFCdata=flattenFC(FCtrain);
            [scores_ACSC] = edgeLearnACSC(nmEdges, size(FCtrain,1), nmScans, flatFCdata);   % ACSC
            scores_ACSC(isnan(scores_ACSC))=0;
            [scores_ACSC_sorted,I1]=sort(scores_ACSC,'descend');
            scores_ACSC_all(:,r,i)=scores_ACSC;
            
            FeaturesNum=[1000:3000:19000];
            for j=1:length(FeaturesNum)          
                % construct train data and test data
                Traindata=[]; Testdata=[];
                for iscan=1:nmScans
                    fc=cell2mat(flatFCdata(iscan));
                    Traindata=[Traindata;fc(I1(1:FeaturesNum(j)),:)'];
                end
                %%% choose features by ACSC
                flatFCdata_test=flattenFC(FCtest);
                for iscan=1:nmScans
                    fc=cell2mat(flatFCdata_test(iscan));
                    Testdata=[Testdata;fc(I1(1:FeaturesNum(j)),:)'];
                end
                % SVM training for multi-class
                svmModel = fitcecoc(Traindata,Ytrain);
                % SVM prediction
                Ypred = predict(svmModel, Testdata);
                % Accuracy and F1 score
                acc_fold_acsc(i,j) = mean(Ypred == Ytest);
                f1_fold_acsc(i,j)=Macro_F1_score(double(Ypred),double(Ytest));
                
            end
            
            % using full fc
            Traindata=[]; Testdata=[];   % using full fc
            for iscan=1:nmScans
                fc=cell2mat(flatFCdata(iscan));
                Traindata=[Traindata;fc'];
            end
            for iscan=1:nmScans
                fc=cell2mat(flatFCdata_test(iscan));
                Testdata=[Testdata;fc'];
            end
            % SVM training for multi-class
            svmModel = fitcecoc(Traindata,Ytrain);
            % SVM prediction
            Ypred = predict(svmModel, Testdata);
            % Accuracy and F1 score
            acc_fold_full(i) = mean(Ypred == Ytest);
            f1_fold_full(i)=Macro_F1_score(double(Ypred),double(Ytest));
        end
        all_accuracies_acsc(ispar,r,:,:)=acc_fold_acsc;
        all_f1_scores_acsc(ispar,r,:,:)=f1_fold_acsc;
        all_accuracies_full(ispar,r,:)=acc_fold_full;
        all_f1_scores_full(ispar,r,:)=f1_fold_full;
    end
    ACSC(ispar-3,:)=mean(mean(scores_ACSC_all,3),2);
    
end
save ClassifyResluts_par.mat ACSC all_accuracies_acsc all_f1_scores_acsc all_accuracies_full all_f1_scores_full









