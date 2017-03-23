function [dist,minX, minY, minZ] = getDistToClosestVoxel(src_vox, target, targetVoxelSize, targetIdx, X, Y, Z)
    dist  = 1e38; 
    minX = -1;
    minY = -1;
    minZ = -1;
    
    for t = 1:size(targetIdx)
        curr_vox = [X(t), Y(t), Z(t)];
        diffVox = (src_vox-curr_vox).*targetVoxelSize;
        currDist = sum(diffVox.*diffVox);
        if currDist < dist
            dist = currDist;
            minX = X(t);
            minY = Y(t);
            minZ = Z(t);
            if dist == 0
                %disp('Found 0')
                break
            end
        end 
    end
    
    
    %{
    for x = 1:s(1)
        for y = 1:s(2)
            for z = 1:s(3)
                curr_vox = [x, y, z];
                curr_vox_intensity = target(x,y,z);
                
                if (curr_vox_intensity>0)
                    %disp(src_vox)
                    %disp(curr_vox)
                    diffVox = (src_vox-curr_vox).*targetVoxelSize;
                    %disp(diffVox)
                    diffVox = diffVox.*diffVox;
                    currDist = sqrt(sum(diffVox));
                    %disp(currDist)
                    if currDist < dist
                        dist = currDist;
                        minX = x;
                        minY = y;
                        minZ = z;
                    end 
                end
            end
        end
    end
    %}
end
