%% FC analysis --- ANOVA  - age, gender ,clip effect on fc
clear;clc
result_path='F:\Projects\20230919_MovieGraphtheory\result\FCA';
mkdir(result_path)
cd(result_path)
info_path='F:\Projects\20230919_MovieGraphtheory\result\SubInfo.xlsx';
[~,~,raw]=xlsread(info_path);
raw(1,:)=[];
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
            info(isub,1)=raw(sub_label,2);   % gender
            info(isub,2)=raw(sub_label,3);   % age
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
for ispar=4:6
    parfor iroi=1:200
        iroi
        for jroi=1:200
            FC=squeeze(FC_all(:,:,ispar,iroi,jroi));
            g1=repmat(info(:,1),1,size(FC,2));
            g2=repmat(info(:,2),1,size(FC,2));
            g3=repmat([1:15],size(FC,1),1);
            group={g1(:),g2(:),g3(:)};
            %%% gender and age difference   -- anovan
            [p,tbl,stats]=anovan(FC(:),group,'model','full','display','off',...
                'varnames', {'Gender', 'AgeGroup', 'Clip'});
            
            gender_p(iroi,jroi)=tbl{2,7};   % main effect -gender
            gender_F(iroi,jroi)=tbl{2,6};
            c=multcompare(stats,'Dimension',[1],'Display','off');
            gender_compare(iroi,jroi)=c(1,4);   % F v.s. M - pos means F greater
            
            age_p(iroi,jroi)=tbl{3,7};   % main effect -age
            age_F(iroi,jroi)=tbl{3,6};
            c=multcompare(stats,'Dimension',[2],'Display','off');
            cc=c(:,4);   % age group compare
            if cc(1)<0 && cc(2)>0 && cc(3) >0
            age_compare(iroi,jroi)=1;   % means elders are greater
            elseif cc(1)>0 && cc(2)<0 && cc(3) <0
                age_compare(iroi,jroi)=-1;  % means elders are smaller
            else
                age_compare(iroi,jroi)=0;
            end
            
            
            clip_p(iroi,jroi)=tbl{4,7};   % main effect -clip
            clip_F(iroi,jroi)=tbl{4,6};
            c=multcompare(stats,'Dimension',[3],'Display','off');
            cc=find(c(:,6)<0.05);
            if length(cc)>=62 
                clip_compare(iroi,jroi)=1;  % content-specific
            else
                clip_compare(iroi,jroi)=0;
            end

            inter_GenderAndAge_p(iroi,jroi)=tbl{5,7};   % INTER
            inter_GenderAndAge_F(iroi,jroi)=tbl{5,6};
            inter_GenderAndClip_p(iroi,jroi)=tbl{6,7};   % INTER
            inter_GenderAndClip_F(iroi,jroi)=tbl{6,6};
            inter_AgeAndClip_p(iroi,jroi)=tbl{7,7};   % INTER
            inter_AgeAndClip_F(iroi,jroi)=tbl{7,6};
        end
    end
    %%% fdr
    [Anova_result.gender_p_fdr(:,:,ispar-3),Anova_result.gender_F_fdr(:,:,ispar-3)]=y_FDR(gender_p,gender_F);
    [~,Anova_result.gender_compare_fdr(:,:,ispar-3)]=y_FDR(gender_p,gender_compare);
    
    [Anova_result.age_p_fdr(:,:,ispar-3),Anova_result.age_F_fdr(:,:,ispar-3)]=y_FDR(age_p,age_F);
    [~,Anova_result.age_compare_fdr(:,:,ispar-3)]=y_FDR(age_p,age_compare);
    
    [Anova_result.clip_p_fdr(:,:,ispar-3),Anova_result.clip_F_fdr(:,:,ispar-3)]=y_FDR(clip_p,clip_F);
    [~,Anova_result.clip_compare_fdr(:,:,ispar-3)]=y_FDR(clip_p,clip_compare);
    
    [Anova_result.inter_GenderAndAge_p_fdr(:,:,ispar-3),Anova_result.inter_GenderAndAge_F_fdr(:,:,ispar-3)]=y_FDR(inter_GenderAndAge_p,inter_GenderAndAge_F);
    [Anova_result.inter_GenderAndClip_p_fdr(:,:,ispar-3),Anova_result.inter_GenderAndClip_F_fdr(:,:,ispar-3)]=y_FDR(inter_GenderAndClip_p,inter_GenderAndClip_F);
    [Anova_result.inter_AgeAndClip_p_fdr(:,:,ispar-3),Anova_result.inter_AgeAndClip_F_fdr(:,:,ispar-3)]=y_FDR(inter_AgeAndClip_p,inter_AgeAndClip_F);
end
save ANOVA.mat info  Anova_result 

%%  Age, gender distribiution  
clear;clc
info_path='F:\Projects\20230919_MovieGraphtheory\result\FCA\ANOVA.mat';
load(info_path);

for isub=1:length(info)
    %%% gender
    if cell2mat(info(isub,1)) == 'M'
        gender(isub) =1;
    else
        gender(isub)=2;
    end
    %%% age
    if string(info(isub,2)) == "22-25"
        agegroup(isub)=1;
    elseif string(info(isub,2)) == "26-30"
        agegroup(isub)=2;
    elseif string(info(isub,2)) == "31-35"
        agegroup(isub)=3;
    else
        agegroup(isub)=4;
    end
end

a=find(gender==2 & agegroup==3);





