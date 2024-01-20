%% Analysis  --- efficiency  --- Ttest and ANOVA
% find fc --- movie > resting; movie < resting; be different among movie
% clips
clear;clc
result_path='F:\Projects\Movie_graph_theory_2023_9_19\result\FCA';
mkdir(result_path)
cd(result_path)
%%% construct data
sup_path='F:\Projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200';
run_dir=dir(sup_path);
run_dir(1:2)=[];
label=1;
FC_all=zeros(171,19,10,200,200);
for irun=1:length(run_dir)
    tic
    if irun == 2   % it is resting   ---  rfMRI_REST1_7T_PA
        data_path=fullfile(sup_path,run_dir(irun).name,'GraphTheory_sparsity');
        sub_dir=dir(data_path);
        sub_dir(1:2)=[];
        for isub=1:length(sub_dir)
            disp(strcat('loading...RunNum_',num2str(irun),'...SubNum_',num2str(isub)));
            load(fullfile(data_path,sub_dir(isub).name));
            var=fieldnames(cell2mat(sw));
            vv=1:length(var);
            vv(find(vv==14 | vv==16))=[];   % locErand and gErand are null
            x=cell2mat(sw);
            for ispar=1:length(fc_sparsity)
                FC_all(isub,label,ispar,:,:)=cell2mat(fc_sparsity(ispar));
                for ivar=1:length(vv)
                    exp=['GraphTheory.',var{vv(ivar)},'(',num2str(isub),',',num2str(label),',',num2str(ispar),',:)',' = ','x(',num2str(ispar),').',var{vv(ivar)},';'];
                    eval(exp)
                end
            end
            toc
        end
        label=label+1;
    elseif irun > 4
        clip_dir=dir(fullfile(sup_path,run_dir(irun).name,'clips_originvalue'));
        clip_dir(1:2)=[];
        for iclip=1:length(clip_dir)
            disp(strcat('loading...RunNum_',num2str(irun),'...Clip_',num2str(iclip)));
            sub_dir=dir(fullfile(clip_dir(iclip).folder,clip_dir(iclip).name,'GraphTheory_sparsity'));
            sub_dir(1:2)=[];
            for isub=1:length(sub_dir)
                load(fullfile(sub_dir(isub).folder,sub_dir(isub).name));
                x=cell2mat(sw);
                for ispar=1:length(fc_sparsity)
                    FC_all(isub,label,ispar,:,:)=cell2mat(fc_sparsity(ispar));
                    for ivar=1:length(vv)
                        exp=['GraphTheory.',var{vv(ivar)},'(',num2str(isub),',',num2str(label),',',num2str(ispar),',:)',' = ','x(',num2str(ispar),').',var{vv(ivar)},';'];
                        eval(exp)
                    end
                end
            end
            label=label+1;
            toc
        end
    end
end
%%
%%%%%%
%%%%%%  anova
var=fieldnames(GraphTheory);
key_num=[12,80];   % based on binomial distribution (p =0.05)
condition={'rest','a1','a2','a3','a4','a5','b1','b2','b3','b4','c1','c2','c3','c4','c5','d1','d2','d3','d4'};
for ispar=1:length(fc_sparsity)
    %%% FC
    tic
    for iroi=1:200
        disp(strcat('Anova_FC...','ispar=',num2str(ispar),'...iroi=',num2str(iroi)));
        for jroi=1:200
            %%% anova
            x=squeeze(FC_all(:,:,ispar,iroi,jroi));
            [pvalue_anova(iroi,jroi),tbl{iroi,jroi},f_value(iroi,jroi)]=anova_rm1(x,condition,0);  % repeated one-way anova
        end
        toc
    end
    pvalue_anova(isnan(pvalue_anova))=1;
    f_value(isnan(f_value))=0;
    [pvalue_anova_fdr,f_value_fdr]=y_FDR(pvalue_anova,f_value);   % fdr
    % judge
    Judge_0=zeros(200,200);
    Judge_1=zeros(200,200);
    Judge_2=zeros(200,200);
    A=diag(repmat([1],19,1));
    tri_label=find(A==0);   % align the tbl to a diag
    sig_label=find(pvalue_anova_fdr<0.05);
    %%% Judge 0  -- the anova result
    Judge_0(sig_label)=1;
    for ilabel=1:length(sig_label)
        tbl_now=tbl{sig_label(ilabel)};
        %%% Judge 1 --- rest v.s clips
        anovaP_now=ones(19,19);
        anovaP_now(tri_label)=tbl_now.pValue;   % current p value
        anovaD_now=zeros(19,19);
        anovaD_now(tri_label)=tbl_now.Difference;   % current deffenrece
        
        kk=sum(anovaP_now(1,2:end)<0.05);  % the significant number -- rest v.s movie clips
        if kk > key_num(1)
            if sum(anovaD_now(1,2:end)<0) > key_num(1)
                Judge_1(sig_label(ilabel))=-1;   % the rest < clips
            elseif sum(anovaD_now(1,2:end)>0) > key_num(1)
                Judge_1(sig_label(ilabel))=1;   % the test > clips
            end
        else
            Judge_1(sig_label(ilabel))=0;
        end
        %%% Judge 2  ---  among clips
        anovaP_now_clips=triu(anovaP_now(2:end,2:end),1);
        kk=sum(sum(anovaP_now_clips<0.05 & anovaP_now_clips >0));  % the significant number -- among 18 clips
        if kk > key_num(2)
            Judge_2(sig_label(ilabel))=1;
        else
            Judge_2(sig_label(ilabel))=0;
        end
    end
    exp=['Result_FC_anova','(:,:,',num2str(ispar),')',' = ','Judge_0',';'];
    eval(exp)
    exp=['Result_FC_restVSclips','(:,:,',num2str(ispar),')',' = ','Judge_1',';'];
    eval(exp)
    exp=['Result_FC_Clips','(:,:,',num2str(ispar),')',' = ','Judge_2',';'];
    eval(exp)
    clear pvalue_anova tbl f_value Judge_0 Judge_1 Judge_2
    
    %%% Graph theory
    for ivar=1:10
        disp(strcat('Anova_Graph...','ispar=',num2str(ispar),'...ivar=',num2str(ivar)));
        exp=['data',' = ','GraphTheory.',var{vv(ivar)},';'];
        eval(exp)
        for iroi=1:200
            %%% anova
            x=squeeze(data(:,:,ispar,iroi));
            [pvalue_anova(iroi),tbl{iroi},f_value(iroi)]=anova_rm1(x,condition,0);  % repeated one-way anova
        end
        pvalue_anova(isnan(pvalue_anova))=1;
        f_value(isnan(f_value))=0;
        [pvalue_anova_fdr,f_value_fdr]=y_FDR(pvalue_anova,f_value);   % fdr
        % judge
        Judge_0=zeros(200,1);
        Judge_1=zeros(200,1);
        Judge_2=zeros(200,1);
        A=diag(repmat([1],19,1));
        tri_label=find(A==0);   % align the tbl to a diag
        sig_label=find(pvalue_anova_fdr<0.05);
        %%% Judge 0  -- the anova result
        Judge_0(sig_label)=1;
        for ilabel=1:length(sig_label)
            tbl_now=tbl{sig_label(ilabel)};
            %%% Judge 1 --- rest v.s clips
            anovaP_now=ones(19,19);
            anovaP_now(tri_label)=tbl_now.pValue;   % current p value
            anovaD_now=zeros(19,19);
            anovaD_now(tri_label)=tbl_now.Difference;   % current deffenrece

            kk=sum(anovaP_now(1,2:end)<0.05);  % the significant number -- rest v.s movie clips
            if kk > key_num(1)
                if sum(anovaD_now(1,2:end)<0) > key_num(1)
                    Judge_1(sig_label(ilabel))=-1;   % the rest < clips
                elseif sum(anovaD_now(1,2:end)>0) > key_num(1)
                    Judge_1(sig_label(ilabel))=1;   % the test > clips
                end
            else
                Judge_1(sig_label(ilabel))=0;
            end
            %%% Judge 2  ---  among clips
            anovaP_now_clips=triu(anovaP_now(2:end,2:end),1);
            kk=sum(sum(anovaP_now_clips<0.05 & anovaP_now_clips >0));  % the significant number -- among 18 clips
            if kk > key_num(2)
                Judge_2(sig_label(ilabel))=1;
            else
                Judge_2(sig_label(ilabel))=0;
            end
        end
        exp=['Result_Graph_anova.',var{ivar},'(',num2str(ispar),',:)',' = ','Judge_0',';'];
        eval(exp)
        exp=['Result_Graph_restVSclips.',var{ivar},'(',num2str(ispar),',:)',' = ','Judge_1',';'];
        eval(exp)
        exp=['Result_Graph_Clips.',var{ivar},'(',num2str(ispar),',:)',' = ','Judge_2',';'];
        eval(exp)
        clear pvalue_anova tbl f_value Judge_0 Judge_1 Judge_2
        toc
    end
    
    for ivar=11:21
        disp(strcat('Anova_Graph...','ispar=',num2str(ispar),'...ivar=',num2str(ivar)));
        exp=['data',' = ','GraphTheory.',var{ivar},';'];
        eval(exp)
        %%% anova
        x=squeeze(data(:,:,ispar));
        [pvalue_anova,tbl,f_value]=anova_rm1(x,condition,0);  % repeated one-way anova
        if pvalue_anova<0.05
            Judge_0=1;
        end
        % judge
        A=diag(repmat([1],19,1));
        tri_label=find(A==0);   % align the tbl to a diag
        tbl_now=tbl;
        %%% Judge 1 --- rest v.s clips
        anovaP_now=ones(19,19);
        anovaP_now(tri_label)=tbl_now.pValue;   % current p value
        anovaD_now=zeros(19,19);
        anovaD_now(tri_label)=tbl_now.Difference;   % current deffenrece

        kk=sum(anovaP_now(1,2:end)<0.05);  % the significant number -- rest v.s movie clips
        if kk > key_num(1)
            if sum(anovaD_now(1,2:end)<0) > key_num(1)
                Judge_1=-1;   % the rest < clips
            elseif sum(anovaD_now(1,2:end)>0) > key_num(1)
                Judge_1=1;   % the test > clips
            else
               Judge_1=0; 
            end
        else
            Judge_1=0;
        end
        
        %%% Judge 2  ---  among clips
        anovaP_now_clips=triu(anovaP_now(2:end,2:end),1);
        kk=sum(sum(anovaP_now_clips<0.05 & anovaP_now_clips >0));  % the significant number -- among 18 clips
        if kk > key_num(2)
            Judge_2=1;
        else
            Judge_2=0;
        end
        Judge_0=(Judge_1 | Judge_2);
        exp=['Result_Graph_anova.',var{ivar},'(',num2str(ispar),')',' = ','Judge_0',';'];
        eval(exp)
        exp=['Result_Graph_restVSclips.',var{ivar},'(',num2str(ispar),')',' = ','Judge_1',';'];
        eval(exp)
        exp=['Result_Graph_Clips.',var{ivar},'(',num2str(ispar),')',' = ','Judge_2',';'];
        eval(exp)
        clear pvalue_anova tbl f_value Judge_0 Judge_1 Judge_2
        toc
    end
end
save FCA_sparsity_fdr.mat Result_FC_anova Result_FC_restVSclips Result_FC_Clips Result_Graph_anova Result_Graph_restVSclips Result_Graph_Clips


%%  FCA   - gender difference
clc;clear
load('F:\Projects\Movie_graph_theory_2023_9_19\result\FCA\FCA_sparsity.mat','FC_all','GraphTheory');
[~,~,raw]=xlsread('F:\Projects\Movie_graph_theory_2023_9_19\result\SubInfo.xlsx');

%%   classify
clear;clc
result_path='F:\Projects\Movie_graph_theory_2023_9_19\result\FCA';
mkdir(result_path)
cd(result_path)
%%% construct data
sup_path='F:\Projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200';
run_dir=dir(sup_path);
run_dir(1:2)=[];
label=1;
FC_all=zeros(171,19,10,200,200);
for irun=1:length(run_dir)
    tic
    if irun == 1   % it is resting   ---  rfMRI_REST1_7T_PA
        data_path=fullfile(sup_path,run_dir(irun).name,'GraphTheory_sparsity');
        sub_dir=dir(data_path);
        sub_dir(1:2)=[];
        for isub=1:length(sub_dir)
            disp(strcat('loading...RunNum_',num2str(irun),'...SubNum_',num2str(isub)));
            load(fullfile(data_path,sub_dir(isub).name));
            var=fieldnames(cell2mat(sw));
            vv=1:length(var);
            vv(find(vv==14 | vv==16))=[];   % locErand and gErand are null
            x=cell2mat(sw);
            for ispar=1:length(fc_sparsity)
                FC_all(isub,label,ispar,:,:)=cell2mat(fc_sparsity(ispar));
                for ivar=1:length(vv)
                    exp=['GraphTheory.',var{vv(ivar)},'(',num2str(isub),',',num2str(label),',',num2str(ispar),',:)',' = ','x(',num2str(ispar),').',var{vv(ivar)},';'];
                    eval(exp)
                end
            end
            toc
        end
        label=label+1;
    elseif irun > 4
        clip_dir=dir(fullfile(sup_path,run_dir(irun).name,'clips_originvalue'));
        clip_dir(1:2)=[];
        for iclip=1:length(clip_dir)
            disp(strcat('loading...RunNum_',num2str(irun),'...Clip_',num2str(iclip)));
            sub_dir=dir(fullfile(clip_dir(iclip).folder,clip_dir(iclip).name,'GraphTheory_sparsity'));
            sub_dir(1:2)=[];
            for isub=1:length(sub_dir)
                load(fullfile(sub_dir(isub).folder,sub_dir(isub).name));
                x=cell2mat(sw);
                for ispar=1:length(fc_sparsity)
                    FC_all(isub,label,ispar,:,:)=cell2mat(fc_sparsity(ispar));
                    for ivar=1:length(vv)
                        exp=['GraphTheory.',var{vv(ivar)},'(',num2str(isub),',',num2str(label),',',num2str(ispar),',:)',' = ','x(',num2str(ispar),').',var{vv(ivar)},';'];
                        eval(exp)
                    end
                end
            end
            label=label+1;
            toc
        end
    end
end
%%%%%%%  Classify
load('F:\Projects\Movie_graph_theory_2023_9_19\result\FCA\FCA_sparsity');
N=200;   % number of rois
nmEdges=(N*(N-1))/2;   % number of edges
nmSbj=171;   % number of subjects
nmScans=19;   % number of scans
k_fold=3;
rounds=1;
method='svr';
var=fieldnames(GraphTheory);
graphIndex=[1,2,3,5,7,9];
%%% using FC
for ispar=1:10
    flatFCdata=cell(1,nmScans);
    %%% choose features by ACSC
    for iscan=1:nmScans
        for isub=1:nmSbj   % make the N*N FC matrix to nmEdges*1 vector
            FCpnts(:,isub)  = flattenMatrixToVectorFC( squeeze(FC_all(isub,iscan,ispar,:,:)), N, 3 );
        end
        flatFCdata(iscan)={FCpnts};
    end
    [scores_ACSC] = edgeLearnACSC(nmEdges, nmSbj, nmScans, flatFCdata);   % ACSC
    scores_ACSC(isnan(scores_ACSC))=0;
    %     [scores_RS] = edgeLearnRS(nmEdges, nmSbj, nmScans, flatFCdata);   % RS
    %     scores_RS(isnan(scores_RS))=0;
    [scores_ACSC_sorted,I1]=sort(scores_ACSC,'descend');
    
    FeaturesNum=[200:100:1000,1500:500:10000];    
    for j=1:length(FeaturesNum)
        timeNow=clock;
        disp(strcat(num2str(timeNow(1)),'-',num2str(timeNow(2)),'-',num2str(timeNow(3)),';;;',...
            num2str(timeNow(4)),':',num2str(timeNow(5)),':',num2str(timeNow(6))));
        data=[];
        for iscan=1:nmScans
            fc=cell2mat(flatFCdata(iscan));
            data=[data;fc(I1(1:FeaturesNum(j)),:)'];
        end
        Y=categorical(reshape(repmat([1:19],nmSbj,1),[],1));
        ACC_fc_acsc1(j)=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - all scans
        
        data=[];
        for iscan=2:nmScans
            fc=cell2mat(flatFCdata(iscan));
            data=[data;fc(I1(1:FeaturesNum(j)),:)'];
        end
        Y=categorical(reshape(repmat([1:18],nmSbj,1),[],1));
        ACC_fc_acsc2(j)=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - all movie clips
        
        data=[];
        for iscan=1:nmScans
            fc=cell2mat(flatFCdata(iscan));
            data=[data;fc(I1(1:FeaturesNum(j)),:)'];
        end
        Y=categorical(reshape(repmat([1,repmat([2],1,18)],nmSbj,1),[],1));
        ACC_fc_acsc3(j)=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - rest v.s. movies
    end
    
    Data=reshape(reshape(squeeze(FC_all(:,:,ispar,:,:)),[nmSbj,nmScans,N*N]),[],N*N);  % fc of all sub/ scan/edges
    %%% using sig edges
    label=find(Result_FC_anova(:,:,ispar)~=0);
    data=Data(:,label);
    Y=categorical(reshape(repmat([1:19],nmSbj,1),[],1));
    ACC_fc_sig1=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - all scans
    
    label=find(Result_FC_Clips(:,:,ispar)~=0);
    data=Data(:,label);
    Y=categorical(reshape(repmat([2:19],nmSbj,1),[],1));
    ACC_fc_sig2=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - all scans
    
    label=find(Result_FC_restVSclips(:,:,ispar)~=0);
    data=Data(:,label);
    Y=categorical(reshape(repmat([1,repmat([2],1,18)],nmSbj,1),[],1));
    ACC_fc_sig3=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - all scans
end
save ClassifyResluts_fc.mat ACC_fc_acsc1 ACC_fc_acsc2 ACC_fc_acsc3 ...
                            ACC_fc_sig1 ACC_fc_sig2 ACC_fc_sig3

%%% using Graph theory index
for ispar=1:10
    for iInd=1:length(graphIndex)
        disp(strcat('sparNum=',num2str(ispar),';;;graphIndex =',var{graphIndex(iInd)}))
        timeNow=clock;
        disp(strcat(num2str(timeNow(1)),'-',num2str(timeNow(2)),'-',num2str(timeNow(3)),';;;',...
            num2str(timeNow(4)),':',num2str(timeNow(5)),':',num2str(timeNow(6))));
        exp=['Data =','GraphTheory.',var{graphIndex(iInd)},';'];
        eval(exp);
        for j =1:2    % - all scans
            switch j
                case 1
                    data=reshape(squeeze(Data(:,:,ispar,:)),[],200);
                case 2
                    exp=['AnovaResult=','Result_Graph_anova.',var{graphIndex(iInd)},';'];
                    eval(exp);
                    label=find(AnovaResult(ispar,:)~=0);
                    data=reshape(squeeze(Data(:,:,ispar,label)),[],length(label));
            end
            Y=categorical(reshape(repmat([1:19],nmSbj,1),[],1));
            ACC_graphIndex1(iInd,j,ispar)=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - all scans
        end
        
        for j =1:2    % - movie clips
            switch j
                case 1
                    data=reshape(squeeze(Data(:,2:19,ispar,:)),[],200);
                case 2
                    exp=['AnovaResult=','Result_Graph_Clips.',var{graphIndex(iInd)},';'];
                    eval(exp);
                    label=find(AnovaResult(ispar,:)~=0);
                    data=reshape(squeeze(Data(:,2:19,ispar,label)),[],length(label));
            end
            Y=categorical(reshape(repmat([2:19],nmSbj,1),[],1));
            ACC_graphIndex2(iInd,j,ispar)=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - movie clips
        end
        
        for j =1:2    % - rest v.s. movies
            switch j
                case 1
                    data=reshape(squeeze(Data(:,:,ispar,:)),[],200);
                case 2
                    exp=['AnovaResult=','Result_Graph_restVSclips.',var{graphIndex(iInd)},';'];
                    eval(exp);
                    label=find(AnovaResult(ispar,:)~=0);
                    data=reshape(squeeze(Data(:,:,ispar,label)),[],length(label));
            end
            Y=categorical(reshape(repmat([1,repmat([2],1,18)],nmSbj,1),[],1));
            ACC_graphIndex3(iInd,j,ispar)=y_MultiClassify(data,Y,k_fold,rounds,method);   % acc - rest v.s. movies
        end
        toc
    end
end
save ClassifyResluts_graphIndex.mat ACC_graphIndex1 ACC_graphIndex2 ACC_graphIndex3

 
