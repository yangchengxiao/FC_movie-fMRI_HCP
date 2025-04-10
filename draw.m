%% FC -- circos map
clear;clc
%%%% construct 200 roi label mat
mask_path='F:\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
% x = cifti_read(mask_path);
x = ft_read_cifti(mask_path,'mapname','array');
Net=["Vis","SomMot","DorsAttn","VentAttn","Limbic","Cont","Default"];
roi_label=[];
roinumOfNet=[];
for inet=1:7
    label{inet}=find(contains(string(x.dlabellabel),Net{inet}));
    roi_label=[roi_label,label{inet}];
    roi_num(inet)=length(label{inet});
    roinumOfNet=[roinumOfNet,1:roi_num(inet)];   % the roi's cumulative number
end
for iroi=1:200   % which net is the roi belong to
    for inet=1:7
        NetCut1=sum(roi_num(1:inet-1))+1;
        NetCut2=sum(roi_num(1:inet))+1;
        if iroi>=NetCut1 && iroi<NetCut2
            roiOfNet(iroi)=inet;
        end
    end
end

%%% reshape matrix according 7 networks
load('F:\Projects\20230919_MovieGraphtheory\result\FCA\ANOVA.mat');
ispar=2;
GenderEffect=Anova_result.gender_F_fdr(:,:,ispar);
GenderEffect2=GenderEffect(roi_label,roi_label);

AgeEffect=Anova_result.age_F_fdr(:,:,ispar);
AgeEffect2=AgeEffect(roi_label,roi_label);

ClipEffect=Anova_result.clip_F_fdr(:,:,ispar);
ClipEffect2=ClipEffect(roi_label,roi_label);

GandA=Anova_result.inter_GenderAndAge_F_fdr(:,:,ispar);
GandA2=GandA(roi_label,roi_label);

GandC=Anova_result.inter_GenderAndClip_F_fdr(:,:,ispar);
GandC2=GandC(roi_label,roi_label);

AandC=Anova_result.inter_AgeAndClip_F_fdr(:,:,ispar);
AandC2=AandC(roi_label,roi_label);

%%%%
A=AgeEffect2;
A=A+diag(inf+zeros(1,length(A)));   % 将对角设置为inf

%% find the fc which are greater in female or olders
clear;clc
%%%% construct 200 roi label mat
mask_path='F:\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
% x = cifti_read(mask_path);
x = ft_read_cifti(mask_path,'mapname','array');
Net=["Vis","SomMot","DorsAttn","VentAttn","Limbic","Cont","Default"];
roi_label=[];
roinumOfNet=[];
for inet=1:7
    label{inet}=find(contains(string(x.dlabellabel),Net{inet}));
    roi_label=[roi_label,label{inet}];
    roi_num(inet)=length(label{inet});
    roinumOfNet=[roinumOfNet,1:roi_num(inet)];   % the roi's cumulative number
end
for iroi=1:200   % which net is the roi belong to
    for inet=1:7
        NetCut1=sum(roi_num(1:inet-1))+1;
        NetCut2=sum(roi_num(1:inet))+1;
        if iroi>=NetCut1 && iroi<NetCut2
            roiOfNet(iroi)=inet;
        end
    end
end
load('F:\Projects\20230919_MovieGraphtheory\result\FCA\ANOVA.mat');
ispar=2;
Gender_compare=Anova_result.gender_compare_fdr(:,:,ispar);
Gender_compare2=Gender_compare(roi_label,roi_label);

Age_compare=Anova_result.age_compare_fdr(:,:,ispar);
Age_compare2=Age_compare(roi_label,roi_label);

k1=1;
k2=1;
for iroi=1:200
    for jroi=1:200
        if Gender_compare2(iroi,jroi)<0    % F v.s. M - pos means F greater
            Gender_comp(k1,1)=roiOfNet(iroi);
            Gender_comp(k1,2)=roiOfNet(jroi);
            k1=k1+1;
        end
        if Age_compare2(iroi,jroi)==-1   % 1 means elders are greater; -1 means elders are smaller
            Age_comp(k2,1)=roiOfNet(iroi);
            Age_comp(k2,2)=roiOfNet(jroi);
            k2=k2+1;
        end
    end
end
for inet=1:7
    label=find(Gender_comp(:,1)==inet & Gender_comp(:,2)==inet);  % intra-network
    exp=['Gender_multcompare.intra_net(',num2str(inet),',1)=',num2str(length(label)),';'];
    eval(exp);
    label=find(Gender_comp(:,1)==inet & Gender_comp(:,2)~=inet);  % inter-network
    exp=['Gender_multcompare.inter_net(',num2str(inet),',1)=',num2str(length(label)),';'];
    eval(exp);
    
    label=find(Age_comp(:,1)==inet & Age_comp(:,2)==inet);  % intra-network
    exp=['Age_multcompare.intra_net(',num2str(inet),',1)=',num2str(length(label)),';'];
    eval(exp);
    label=find(Age_comp(:,1)==inet & Age_comp(:,2)~=inet);  % inter-network
    exp=['Age_multcompare.inter_net(',num2str(inet),',1)=',num2str(length(label)),';'];
    eval(exp);
end

%% find the fc which is sig clip effect in anova
%%% Find the percentage of all edges where the effect of clip on FC is significant.
clear;clc
%%%% construct 200 roi label mat
mask_path='F:\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
% x = cifti_read(mask_path);
x = ft_read_cifti(mask_path,'mapname','array');
Net=["Vis","SomMot","DorsAttn","VentAttn","Limbic","Cont","Default"];
roi_label=[];
roinumOfNet=[];
for inet=1:7
    label{inet}=find(contains(string(x.dlabellabel),Net{inet}));
    roi_label=[roi_label,label{inet}];
    roi_num(inet)=length(label{inet});
    roinumOfNet=[roinumOfNet,1:roi_num(inet)];   % the roi's cumulative number
end
for iroi=1:200   % which net is the roi belong to
    for inet=1:7
        NetCut1=sum(roi_num(1:inet-1))+1;
        NetCut2=sum(roi_num(1:inet))+1;
        if iroi>=NetCut1 && iroi<NetCut2
            roiOfNet(iroi)=inet;
        end
    end
end
load('F:\Projects\20230919_MovieGraphtheory\result\FCA\ANOVA.mat');
ispar=2;
Clip_effect=Anova_result.clip_F_fdr(:,:,ispar);
Clip_effect2=Clip_effect(roi_label,roi_label);

Clip_inter_gender=Anova_result.inter_GenderAndClip_F_fdr(:,:,ispar);
Clip_inter_gender2=Clip_inter_gender(roi_label,roi_label);

Clip_inter_age=Anova_result.inter_AgeAndClip_F_fdr(:,:,ispar);
Clip_inter_age2=Clip_inter_age(roi_label,roi_label);

k1=1;k2=1;k3=1; k4=1;
all_fc=ones(200,200);
all_fc(all_fc==eye(200))=0;
for iroi=1:200
    for jroi=1:200
        if Clip_effect2(iroi,jroi)~=0
            Clip(k1,1)=roiOfNet(iroi);
            Clip(k1,2)=roiOfNet(jroi);
            k1=k1+1;
        end
        if Clip_inter_gender2(iroi,jroi)~=0
            ClipGender(k2,1)=roiOfNet(iroi);
            ClipGender(k2,2)=roiOfNet(jroi);
            k2=k2+1;
        end
        if Clip_inter_age2(iroi,jroi)~=0
            ClipAge(k3,1)=roiOfNet(iroi);
            ClipAge(k3,2)=roiOfNet(jroi);
            k3=k3+1;
        end
        if all_fc(iroi,jroi)~=0
            ALL_edge(k4,1)=roiOfNet(iroi);
            ALL_edge(k4,2)=roiOfNet(jroi);
            k4=k4+1;
        end
    end
end
for inet=1:7
    label1=find(Clip(:,1)==inet & Clip(:,2)==inet);  % intra-network
    label_all1=find(ALL_edge(:,1)==inet & ALL_edge(:,2)==inet);
    label2=find(Clip(:,1)==inet & Clip(:,2)~=inet);  % inter-network
    label_all2=find(ALL_edge(:,1)==inet & ALL_edge(:,2)~=inet);
    Origin_mat_clip(inet,1)=length(label1);   % intra  - sig
    Origin_mat_clip(inet,2)=length(label_all1)-length(label1);   % no-sig
    Origin_mat_clip(inet,3)=length(label2);   % inter  - sig
    Origin_mat_clip(inet,4)=length(label_all2)-length(label2);   % no-sig
    
    label1=find(ClipGender(:,1)==inet & ClipGender(:,2)==inet);  % intra-network
    label2=find(ClipGender(:,1)==inet & ClipGender(:,2)~=inet);  % inter-network
    Origin_mat_ClipGender(inet,1)=length(label1);   % intra  - sig
    Origin_mat_ClipGender(inet,2)=length(label_all1)-length(label1);   % no-sig
    Origin_mat_ClipGender(inet,3)=length(label2);   % inter  - sig
    Origin_mat_ClipGender(inet,4)=length(label_all2)-length(label2);   % no-sig
    
    label1=find(ClipAge(:,1)==inet & ClipAge(:,2)==inet);  % intra-network
    label2=find(ClipAge(:,1)==inet & ClipAge(:,2)~=inet);  % inter-network
    Origin_mat_ClipAge(inet,1)=length(label1);   % intra  - sig
    Origin_mat_ClipAge(inet,2)=length(label_all1)-length(label1);   % no-sig
    Origin_mat_ClipAge(inet,3)=length(label2);   % inter  - sig
    Origin_mat_ClipAge(inet,4)=length(label_all2)-length(label2);   % no-sig
end
%%  find the top 1000 ACSC distribution 
%%% and calculate the mean score on the network level
clear;clc
%%%% construct 200 roi label mat
mask_path='F:\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
% x = cifti_read(mask_path);
x = ft_read_cifti(mask_path,'mapname','array');
Net=["Vis","SomMot","DorsAttn","VentAttn","Limbic","Cont","Default"];
roi_label=[];
roinumOfNet=[];
for inet=1:7
    label{inet}=find(contains(string(x.dlabellabel),Net{inet}));
    roi_label=[roi_label,label{inet}];
    roi_num(inet)=length(label{inet});
    roinumOfNet=[roinumOfNet,1:roi_num(inet)];   % the roi's cumulative number
end
for iroi=1:200   % which net is the roi belong to
    for inet=1:7
        NetCut1=sum(roi_num(1:inet-1))+1;
        NetCut2=sum(roi_num(1:inet))+1;
        if iroi>=NetCut1 && iroi<NetCut2
            roiOfNet(iroi)=inet;
        end
    end
end
load('F:\Projects\20230919_MovieGraphtheory\result\FCA\ClassifyResluts_par.mat');
scores_ACSC=ACSC(2,:);
[scores_ACSC_sorted,I1]=sort(scores_ACSC,'descend');

aa=ones(200,200);
tMask = triu(true(size(aa)), 1);
k=1;
for jroi=1:200
    for iroi=1:200
        if tMask(iroi,jroi)~=0
            roiLabel(k,1)=iroi;
            roiLabel(k,2)=jroi;
            k=k+1;
        end
    end
end
label_top=I1(1:1000);
roiLabel_top=roiLabel(label_top,:);
scores=scores_ACSC_sorted(label_top);
for ilabel=1:length(roiLabel_top)
    Net1=roiOfNet(find(roi_label==roiLabel_top(ilabel,1)));
    Net2=roiOfNet(find(roi_label==roiLabel_top(ilabel,2)));
    RoiNet_top(ilabel,1)=Net1;
    RoiNet_top(ilabel,2)=Net2;
end

for inet=1:7
    label=find(RoiNet_top(:,1)==inet & RoiNet_top(:,2)==inet);  % intra-network
    exp=['ACSC_top.intra_net_',char(Net(inet)),'=',num2str(length(label)),';'];
    eval(exp);
    Orogin_ACSCscore(inet,1)=length(label); 
    label=find(RoiNet_top(:,1)==inet & RoiNet_top(:,2)~=inet);  % inter-network
    exp=['ACSC_top.inter_net_',char(Net(inet)),'=',num2str(length(label)),';'];
    eval(exp);
    Orogin_ACSCscore(inet,2)=length(label);
end





