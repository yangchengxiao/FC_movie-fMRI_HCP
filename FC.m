%% Calculate FC --- the whole time course  2023.9.22
clear;clc;close all
run_path='E:\DataBase\HCP_movie_7T\Preproc_2mm32k_FIX-Denoised';
run_dir=dir(run_path);
run_dir(1:2)=[];

mask_path='D:\ycx\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
out_path='D:\ycx\projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200_GSR';

for irun=1:8
    irun
    sub_path=fullfile(run_path,run_dir(irun).name);
    result_path=fullfile(out_path,run_dir(irun).name,'whole_originvalue');
    mkdir(result_path)
    y_staticFC(sub_path,mask_path,result_path)
end

% Calculate FC  --- divide movie run to clips 2023.9.22
clear;clc;close all
run_path='E:\DataBase\HCP_movie_7T\Preproc_2mm32k_FIX-Denoised';
run_dir=dir(run_path);
run_dir(1:2)=[];

mask_path='D:\ycx\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
out_path='D:\ycx\projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200_GSR';

clips_label1=[[21,274];[285,516];[527,724];[735,808];[819,911]];  % The start and end time_point of the clips of movie1
clips_label2=[[21,257];[268,536];[547,805];[816,898]];    % The start and end time_point of the clips of movie2
clips_label3=[[21,211];[222,416];[427,640];[651,802];[813,905]];    % The start and end time_point of the clips of movie3
clips_label4=[[21,263];[274,513];[524,788];[799,891]];    % The start and end time_point of the clips of movie4
clips_num=[5,4,5,4];   % the clips number of each movie

for irun=5:8
    sub_path=fullfile(run_path,run_dir(irun).name);
    sub_dir=dir(sub_path);
    sub_dir(1:2)=[];
    for isub=1:length(sub_dir)
        tic
        disp(strcat('irun = ',num2str(irun-4),'...isub = ',num2str(isub)));
        nii_dir=dir(fullfile(sub_path,sub_dir(isub).name,'*MSMAll_hp2000_clean.dtseries.nii'));
        if length(nii_dir)==0
            continue   % if there is no such a nii file, skip the subject
        end
        BrainData=cifti_read(fullfile(nii_dir.folder,nii_dir.name));
        BrainMask=cifti_read(mask_path);
        tp=size(BrainData.cdata,2);
        vert_num_left=BrainData.diminfo{1, 1}.models{1, 1}.numvert;
        vert_num_right=BrainData.diminfo{1, 1}.models{1, 2}.numvert;
        % extract surf_data from brain_data
        surf_data=zeros(length(BrainMask.cdata),tp);
        k_left=BrainData.diminfo{1, 1}.models{1, 1}.vertlist+1;
        surf_data(k_left,:)=BrainData.cdata(BrainData.diminfo{1, 1}.models{1, 1}.start:BrainData.diminfo{1, 1}.models{1, 2}.start-1,:);
        k_right=BrainData.diminfo{1, 1}.models{1, 2}.vertlist+1+vert_num_left;
        surf_data(k_right,:)=BrainData.cdata(BrainData.diminfo{1, 1}.models{1, 2}.start:BrainData.diminfo{1, 1}.models{1, 3}.start-1,:);
        for iclip=1:clips_num(irun-4)
            expr=['clips_now','= ','clips_label',num2str(irun-4),'(iclip,:);'];
            eval(expr);   % The start and end time_point of current clip
            surf_data_clip=surf_data(:,clips_now(1):clips_now(2));
            tp_clip=size(surf_data_clip,2);
            % calculate ROIS
            mask_label=unique(BrainMask.cdata);
            rois=zeros(length(mask_label)-1,tp_clip);
            rois_gsr=rois;
            gs=mean(surf_data_clip(find(BrainMask.cdata ~= 0),:),1)';  % gs
            for iroi=1:length(mask_label)-1   % rois
                mask_location=find(BrainMask.cdata==iroi);
                rois(iroi,:)=mean(surf_data_clip(mask_location,:));
                [~,~,rois_gsr(iroi,:)]=regress(gs,[ones(size(gs,1),1) rois(iroi,:)']);
            end
            [R,p]=corr(rois');
            R=(R+R')/2;
            R(isnan(R))=0;
            R=R-diag(diag(R));
            %    [~,R]=y_FDR(p,R);
            Z=0.5.*log((1+R)./(1-R));
            
            %%% GSR
            [R_gsr,p]=corr(rois_gsr');
            R_gsr=(R_gsr+R_gsr')/2;
            R_gsr(isnan(R_gsr))=0;
            R_gsr=R_gsr-diag(diag(R_gsr));
            %    [~,R_gsr]=y_FDR_gsr(p,R_gsr);
            Z_gsr=0.5.*log((1+R_gsr)./(1-R_gsr));
            
            out_putR=fullfile(out_path,run_dir(irun).name,'clips_originvalue',['clip',num2str(iclip)],'StaticFC');
            mkdir(out_putR);
            FCFileR=fullfile(out_putR,sub_dir(isub).name);
            save(FCFileR,'R','Z','R_gsr','Z_gsr');
        end
        toc
    end
end

%% graph theory ---    REST 4 RUNS  --- whole time course
clear;clc;close all
run_path='D:\ycx\projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200_GSR';
run_dir=dir(run_path);
run_dir(1:2)=[];

label_exclude=xlsread('D:\ycx\projects\Movie_graph_theory_2023_9_19\result\delete sub.xlsx');  % load label of the exclude

rand_n=1000;
for irun=1:4
    result_path=fullfile(run_path,run_dir(irun).name,'GraphTheory_sparsity_gsr');
    mkdir(result_path);
    cd(result_path)
    
    fc_dir=dir(fullfile(run_path,run_dir(irun).name,'whole_originvalue','StaticFC'));
    fc_dir(1:2)=[];
    %%% delete subjects
    kk=1;
    label=[];
    for isub=1:length(fc_dir)
        a=string(label_exclude);
        if sum(ismember(a,fc_dir(isub).name(1:end-4)))
            label(kk)=isub;
            kk=kk+1;
        end
    end
    fc_dir(label)=[];     % Set the subjects to be excluded as empty
    clear label
    for isub=1:length(fc_dir)
        tic
        disp(strcat('irun =',num2str(irun),'.....isub = ',num2str(isub)));
        load(fullfile(fc_dir(isub).folder,fc_dir(isub).name));
        fc=Z_gsr;
        %         fc_binary=y_binary(fc);   % binarized fc
        sparsitys=0.05:0.05:0.5;
        parfor ispar=1:length(sparsitys)
            [fc_sparsity{ispar}, rthr(ispar)] = y_SparsityMatrix(fc,  sparsitys(ispar));
            [N(ispar), K(ispar), sw{ispar}] = y_GraphTheory(fc_sparsity{ispar}, rand_n);
        end
        fname=fc_dir(isub).name(1:end-4);
        save(fname,'fc','fc_sparsity','rthr','N','K','sw');
        toc
    end
end
%%
%%% graph theory   --- clips
clear;clc;close all
run_path='D:\ycx\projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200';
run_dir=dir(run_path);
run_dir(1:2)=[];

clips_label1=[[21,274];[285,516];[527,724];[735,808];[819,911]];  % The start and end time_point of the clips of movie1
clips_label2=[[21,257];[268,536];[547,805];[816,898]];    % The start and end time_point of the clips of movie2
clips_label3=[[21,211];[222,416];[427,640];[651,802];[813,905]];    % The start and end time_point of the clips of movie3
clips_label4=[[21,263];[274,513];[524,788];[799,891]];    % The start and end time_point of the clips of movie4
clips_num=[5,4,5,4];   % the clips number of each movie

label_exclude=xlsread('D:\ycx\projects\Movie_graph_theory_2023_9_19\result\delete sub.xlsx');  % load label of the exclude

rand_n=1000;
for irun=7:8
    for iclip=1:clips_num(irun-4)
        result_path=fullfile(run_path,run_dir(irun).name,'clips_originvalue',['clip',num2str(iclip)],'GraphTheory_sparsity');
        mkdir(result_path);
        cd(result_path)
        fc_dir=dir(fullfile(run_path,run_dir(irun).name,'clips_originvalue',['clip',num2str(iclip)],'StaticFC'));
        fc_dir(1:2)=[];
        %%% delete subjects
        kk=1;
        label=[];
        for isub=1:length(fc_dir)
            a=string(label_exclude);
            if sum(ismember(a,fc_dir(isub).name(1:end-4)))
                label(kk)=isub;
                kk=kk+1;
            end
        end
        fc_dir(label)=[];     % Set the subjects to be excluded as empty
        clear label
        for isub=1:length(fc_dir)
            tic
            disp(strcat('irun =',num2str(irun-4),'...iclip =',num2str(iclip),'...isub = ',num2str(isub)))
            load(fullfile(fc_dir(isub).folder,fc_dir(isub).name));
            fc=Z;
%             fc_binary=y_binary(fc);   % binarized fc
            
            % network efficiency (local and global efficiency) of a weight FC
            sparsitys=0.05:0.05:0.5;
            parfor ispar=1:length(sparsitys)
                [fc_sparsity{ispar}, rthr(ispar)] = y_SparsityMatrix(fc,  sparsitys(ispar));
                [N(ispar), K(ispar), sw{ispar}] = y_GraphTheory(fc_sparsity{ispar}, rand_n);
            end
            fname=fc_dir(isub).name(1:end-4);
            save(fname,'fc','fc_sparsity','rthr','N','K','sw');
            toc
        end
    end
end


%% Calculate dynamic FC 2024.3.3
clear;clc;close all
run_path='H:\DataBase\HCP_movie_7T\Preproc_2mm32k_FIX-Denoised';
run_dir=dir(run_path);
run_dir(1:2)=[];

mask_path='F:\Mask\HCP\fslr32k_cifti\Schaefer2018_200Parcels_7Networks_order.dlabel.nii';
out_path='H:\projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200_dynamic';
window=50;step=5;
for irun=5:8
    irun
    sub_path=fullfile(run_path,run_dir(irun).name);
    result_path=fullfile(out_path,run_dir(irun).name,['window_',num2str(window),'_step_',num2str(step)]);
    mkdir(result_path)
    
    y_dynamicFC(sub_path,mask_path,result_path,window,step)
end
