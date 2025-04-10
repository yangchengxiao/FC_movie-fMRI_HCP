function   flatFCdata=flattenFC(FC)
[nmSub,nmScans,N,~]=size(FC);
for iscan=1:nmScans
    for isub=1:nmSub   % make the N*N FC matrix to nmEdges*1 vector
        FCpnts(:,isub)  = flattenMatrixToVectorFC( squeeze(FC(isub,iscan,:,:)), N, 3 );
    end
    flatFCdata(iscan)={FCpnts};
end
