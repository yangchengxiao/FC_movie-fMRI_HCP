%%   classify   -- LOOCV
clear;clc
result_path='F:\Projects\20230919_MovieGraphtheory\result\FCA';
mkdir(result_path)
cd(result_path)
%%%%% construct data
sup_path='F:\Projects\20230919_MovieGraphtheory\result\FIX-Denoised_200';
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
%%
%%%%%%%  Classify
N=200;   % number of rois
nmEdges=(N*(N-1))/2;   % number of edges
nmSbj=169;   % number of subjects
nmScans=15;   % number of scans
ytrain=categorical(reshape(repmat([1:15],nmSbj-1,1),[],1));
ytest=categorical(reshape(repmat([1:15],1,1),[],1));
%%% SVM
for ispar=4:6
    h= waitbar(0,'Pls wait......','Name','cross-validation: ...');
    tic
    %%%
    % Perform leave one out cross-validation
    parfor Isub=1:nmSbj
        % Split data into training and test sets
        FCtrain=squeeze(FC_all([1:Isub-1, Isub+1:end],:,ispar,:,:));
        FCtest=squeeze(FC_all(Isub,:,ispar,:,:));
        
        % ACSC feature selection
        flatFCdata=flattenFC(FCtrain);
        [scores_ACSC] = edgeLearnACSC(nmEdges, nmSbj-1, nmScans, flatFCdata);   % ACSC
        scores_ACSC(isnan(scores_ACSC))=0;
        %     [scores_RS] = edgeLearnRS(nmEdges, nmSbj, nmScans, flatFCdata);   % RS
        %     scores_RS(isnan(scores_RS))=0;
        [scores_ACSC_sorted,I1]=sort(scores_ACSC,'descend');
        scores_ACSC_all(Isub,:)=scores_ACSC;
        
        FeaturesNum=[1000:3000:19000];
        predictedLabels_acsc=zeros(nmScans,length(FeaturesNum));
        accuracy_acsc=zeros(length(FeaturesNum),1);
        for j=1:length(FeaturesNum)
            timeNow=clock;
            disp(strcat('sparNum is -',num2str(ispar),'...Isub is',num2str(Isub),'...FeaturesNum is -',num2str(j)));
            disp(strcat(num2str(timeNow(1)),'-',num2str(timeNow(2)),'-',num2str(timeNow(3)),';;;',...
                num2str(timeNow(4)),':',num2str(timeNow(5)),':',num2str(timeNow(6))));
            
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
            svmModel = fitcecoc(Traindata,ytrain);
            % SVM prediction
            predictedLabels_acsc(:,j) = predict(svmModel, Testdata);
            accuracy_acsc(j) = length(find(double(predictedLabels_acsc(:,j))==double(ytest)))/15;
        end
        PredictedLabels_acsc(:,:,Isub)=predictedLabels_acsc;
        Accuracy_acsc(:,Isub)=accuracy_acsc;
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
        svmModel = fitcecoc(Traindata,ytrain);
        % SVM prediction
        predictedLabels_full(:,Isub) = predict(svmModel, Testdata);
        accuracy_full(Isub) = length(find((double(predictedLabels_full(:,Isub)))==double(ytest)))/15;
    end
    
    for j=1:7
        [macro_F1_acsc(j,ispar),micro_F1_acsc(j,ispar)]=compute_f1_scores(repmat([1:15],169,1)', squeeze(PredictedLabels_acsc(:,j,:)));
    end
    [macro_F1_full(ispar),macro_F1_full(ispar)]=compute_f1_scores(repmat([1:15],169,1)', predictedLabels_full);
    
    ACSC(ispar-3,:)=mean(scores_ACSC_all,1);
end
save ClassifyResluts.mat ACSC accuracy_full Accuracy_acsc   ...
    macro_F1_acsc macro_F1_full micro_F1_acsc macro_F1_full

















