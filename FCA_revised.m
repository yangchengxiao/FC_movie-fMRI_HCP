%% FC analysis ---  Linear Mixed-Effects Model   - age, gender ,clip effect on fc
clear;clc
result_path='F:\Projects\20230919_MovieGraphtheory\result\FCA';
mkdir(result_path)
cd(result_path)
% load demo info
info_path='F:\Projects\20230919_MovieGraphtheory\result\SubInfo2.xlsx';
[~,~,raw]=xlsread(info_path);
raw(1,:)=[];
% load head movement data
load('F:\Projects\20230919_MovieGraphtheory\result\HeadMovement.mat');
%%% construct data
sup_path='F:\Projects\20230919_MovieGraphtheory\result\FIX-Denoised_200';
run_dir=dir(sup_path);
run_dir(1:2)=[];
label=1;
FC_all=zeros(169,15,10,200,200);
clips_num=[4,3,4,3];
for irun=5:8
    clip_dir=dir(fullfile(sup_path,run_dir(irun).name,'clips_originvalue'));
    clip_dir(1:2)=[];
    for iclip=1:clips_num(irun-4)
        tic
        disp(strcat('loading...RunNum_',num2str(irun),'...Clip_',num2str(iclip)));
        sub_dir=dir(fullfile(clip_dir(iclip).folder,clip_dir(iclip).name,'GraphTheory_sparsity'));
        sub_dir(1:2)=[];
        for isub=1:length(sub_dir)
            %%% info
            sub_label=find(string(raw(:,1))==sub_dir(isub).name(1:end-4));
            info(isub,1)=cell2mat(raw(sub_label,2));   % gender
            info(isub,2)=cell2mat(raw(sub_label,3));   % age
            load(fullfile(sub_dir(isub).folder,sub_dir(isub).name));
            for ispar=1:length(fc_sparsity)
                FC_all(isub,label,ispar,:,:)=cell2mat(fc_sparsity(ispar));
            end
        end
        label=label+1;
        toc
    end
end
% validation clip  - the average of four movie runs
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
head_movement=zeros(169,15);
head_movement(:,1:4)=squeeze(HeadMovement.RelativeRMS_mean(:,1,1:4));  % movie1-clip1:4
head_movement(:,5:7)=squeeze(HeadMovement.RelativeRMS_mean(:,2,1:3));  % movie2-clip1:3
head_movement(:,8:11)=squeeze(HeadMovement.RelativeRMS_mean(:,3,1:4));  % movie3-clip1:4
head_movement(:,12:14)=squeeze(HeadMovement.RelativeRMS_mean(:,4,1:3));  % movie4-clip1:4
head_movement(:,15)=RelativeRMS_mean_validation;  % mean of movie1-4 - validation clip
HM=reshape(head_movement,[],1);

for ispar=4:6
    parfor iroi=1:200
        timeNow=clock;
        fprintf('-SparNum  %d -- iroi %d ',ispar,iroi)
        disp(strcat(num2str(timeNow(1)),'-',num2str(timeNow(2)),'-',num2str(timeNow(3)),';;;',...
            num2str(timeNow(4)),':',num2str(timeNow(5)),':',num2str(timeNow(6))));
        for jroi=1:200
            if iroi~=jroi
                FC=squeeze(FC_all(:,:,ispar,iroi,jroi));
                Gender=categorical(reshape(repmat(info(:,1),1,size(FC,2)),[],1));  % sex
                AgeGroup=reshape(repmat(info(:,2),1,size(FC,2)),[],1);   % age group
                MovieClip=categorical(reshape(repmat([1:15],size(FC,1),1),[],1));   % % clip number
                SubjectID=categorical(reshape(repmat([1:169],size(FC,2),1)',[],1));
                FC=reshape(FC,[],1);
                %%% gender,age, clip, and hm effects on fc
                tbl = table(FC, AgeGroup, Gender, MovieClip, SubjectID, ...
                    'VariableNames', {'FC', 'Age', 'Gender', 'MovieClip','Subject'});
                
                lme = fitlme(tbl, 'FC ~ Age * Gender * MovieClip + (1|Subject)');
                AncovaResults = anova(lme);
                
                gender_p(iroi,jroi)=AncovaResults{3,5};   % main effect -gender
                gender_F(iroi,jroi)=AncovaResults{3,2};
                gender_compare(iroi,jroi)=lme.Coefficients{3,2};   % F v.s. M - pos means F greater
                gender_R2(iroi,jroi)=Marginal_R2(tbl, 'FC ~ Gender+ (1|Subject)');
                
                age_p(iroi,jroi)=AncovaResults{2,5};   % main effect -age
                age_F(iroi,jroi)=AncovaResults{2,2};
                age_compare(iroi,jroi)=lme.Coefficients{2,2}; % pos means elders are greater
                age_R2(iroi,jroi)=Marginal_R2(tbl, 'FC ~ Age+ (1|Subject)');
                
                clip_p(iroi,jroi)=AncovaResults{4,5};   % main effect -clip
                clip_F(iroi,jroi)=AncovaResults{4,2};
                clip_R2(iroi,jroi)=Marginal_R2(tbl, 'FC ~ MovieClip+ (1|Subject)');
                numClips = 15;
                pvals = ones(numClips, numClips);
                for i = 1:numClips
                    for j = i+1:numClips
                        contrast = zeros(1, length(lme.Coefficients.Name));
                        % Find correct index for MovieClip_i and MovieClip_j
                        clip_i_index = find(strcmp(lme.Coefficients.Name, sprintf('MovieClip_%d', i)));
                        clip_j_index = find(strcmp(lme.Coefficients.Name, sprintf('MovieClip_%d', j)));
                        
                        contrast(clip_i_index) = 1;
                        contrast(clip_j_index) = -1;
                        
                        pvals(i, j) = coefTest(lme, contrast);
                    end
                end
                cc=find(pvals<0.05);
                if length(cc)>=62
                    clip_compare(iroi,jroi)=1;  % content-specific
                else
                    clip_compare(iroi,jroi)=0;
                end
                
                inter_GenderAndAge_p(iroi,jroi)=AncovaResults{5,5};   % INTER
                inter_GenderAndAge_F(iroi,jroi)=AncovaResults{5,2};
                inter_GenderAndAge_R2(iroi,jroi)=Marginal_R2(tbl, 'FC ~ Gender*Age+ (1|Subject)');
                
                inter_GenderAndClip_p(iroi,jroi)=AncovaResults{7,5};   % INTER
                inter_GenderAndClip_F(iroi,jroi)=AncovaResults{7,2};
                inter_GenderAndClip_R2(iroi,jroi)=Marginal_R2(tbl, 'FC ~ Gender*MovieClip + (1|Subject)');
                
                inter_AgeAndClip_p(iroi,jroi)=AncovaResults{6,5};   % INTER
                inter_AgeAndClip_F(iroi,jroi)=AncovaResults{6,2};
                inter_AgeAndClip_R2(iroi,jroi)=Marginal_R2(tbl, 'FC ~ Age*MovieClip + (1|Subject)');
            end
        end
    end
    %%% fdr
    [Ancova_result.gender_p_fdr(:,:,ispar-3),Ancova_result.gender_F_fdr(:,:,ispar-3)]=y_FDR(gender_p,gender_F);
    [~,Ancova_result.gender_compare_fdr(:,:,ispar-3)]=y_FDR(gender_p,gender_compare);
    [~,Ancova_result.gender_R2_fdr(:,:,ispar-3)]=y_FDR(gender_p,gender_R2);
    
    [Ancova_result.age_p_fdr(:,:,ispar-3),Ancova_result.age_F_fdr(:,:,ispar-3)]=y_FDR(age_p,age_F);
    [~,Ancova_result.age_compare_fdr(:,:,ispar-3)]=y_FDR(age_p,age_compare);
    [~,Ancova_result.age_R2_fdr(:,:,ispar-3)]=y_FDR(age_p,age_R2);
    
    [Ancova_result.clip_p_fdr(:,:,ispar-3),Ancova_result.clip_F_fdr(:,:,ispar-3)]=y_FDR(clip_p,clip_F);
    [~,Ancova_result.clip_compare_fdr(:,:,ispar-3)]=y_FDR(clip_p,clip_compare);
    [~,Ancova_result.clip_R2_fdr(:,:,ispar-3)]=y_FDR(clip_p,clip_R2);
    
    [Ancova_result.inter_GenderAndAge_p_fdr(:,:,ispar-3),Ancova_result.inter_GenderAndAge_F_fdr(:,:,ispar-3)]=y_FDR(inter_GenderAndAge_p,inter_GenderAndAge_F);
    [~,Ancova_result.inter_GenderAndAge_R2_fdr(:,:,ispar-3)]=y_FDR(inter_GenderAndAge_p,inter_GenderAndAge_R2);
    
    [Ancova_result.inter_GenderAndClip_p_fdr(:,:,ispar-3),Ancova_result.inter_GenderAndClip_F_fdr(:,:,ispar-3)]=y_FDR(inter_GenderAndClip_p,inter_GenderAndClip_F);
    [~,Ancova_result.inter_GenderAndClip_R2_fdr(:,:,ispar-3)]=y_FDR(inter_GenderAndClip_p,inter_GenderAndClip_R2);
    
    [Ancova_result.inter_AgeAndClip_p_fdr(:,:,ispar-3),Ancova_result.inter_AgeAndClip_F_fdr(:,:,ispar-3)]=y_FDR(inter_AgeAndClip_p,inter_AgeAndClip_F);
    [~,Ancova_result.inter_AgeAndClip_R2_fdr(:,:,ispar-3)]=y_FDR(inter_AgeAndClip_p,inter_AgeAndClip_R2);
end
save ANCOVA.mat info  Ancova_result