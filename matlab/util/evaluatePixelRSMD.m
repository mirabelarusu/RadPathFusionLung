function [rmsd, rmsdStd] = evaluatePixelRSMD(target, src, voxelSize)
    
    meanRmsd = 0;
    countVox = 0;
    
    meanRmsdSq = 0;
    
    s = size(src);
    targetIdx = find(target>0);
    [X,Y,Z] = ind2sub(s,targetIdx);
    
    srcIdx  = find(src>0);
    [X_s,Y_s,Z_s] = ind2sub(s,srcIdx);
    
    for ss = 1:size(srcIdx)
        curr_vox = [X_s(ss), Y_s(ss), Z_s(ss)];
        [sqDist, minX, minY, minZ] = getDistToClosestVoxel(curr_vox, target, voxelSize, targetIdx, X,Y,Z);
        dist = sqrt(sqDist);
        %disp([ss, X_s(ss), Y_s(ss), Z_s(ss), dist, minX, minY, minZ])
        meanRmsd = meanRmsd+dist;
        meanRmsdSq = meanRmsdSq + dist*dist;
    end
    
    if (size(srcIdx,1)>0)
        rmsd = meanRmsd/size(srcIdx,1);
        rmsdStd = sqrt(meanRmsdSq/size(srcIdx,1) - rmsd*rmsd);
    else
        rmsd = 1e38;
        rmsdStd = 1e38;
    end
    
end