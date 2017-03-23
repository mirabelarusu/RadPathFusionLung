function out_vol = createReconstructionTarget(vol, i, target, K, score_forward, rigid_params)

if (score_forward)
    t_slices = 2*K+1;
else
    t_slices = K+1;
end

% create intermediate variables needed in the registration
out_vol = zeros(size(vol,1),size(vol,2), t_slices);
if (i<=K)
    numberSteps = K+2-i;
else
    numberSteps = 1;
end;

% at the K+1th slice use the ex vivo mri if volume exist
if (numel(size(target))>2) % target exist
    out_vol(:,:,K+1) = target(:,:,i)';
end

% fill volumes for registration
added = 0;
for (j = K:-1:numberSteps)
    if (maxall(vol(:,:,i-(K-j+1)))>0 & added < K)
        disp ([num2str(i),' ' , num2str(j), ' ', ...
            num2str(i -(K-j+1)), ' ', num2str(added) ])
        out_vol(:,:,j)    = uint32(vol(:,:,i -(K-j+1)));
        added = added + 1;
    end;
end

if (score_forward==1)
    steps = K;
    if ((i+K)>size(vol,3))
        steps = size(vol,3)-i;
    end
    
    for (j = 1:steps)
        if (maxall(vol(:,:,i+j))>0 & added < 2*K)
            disp ([num2str(i),' ' , num2str(K+1+j), ' ', ...
                num2str(i+j), ' ', num2str(steps) ])
            out_vol(:,:,K+1+j)    = uint32(imdef(rigid_params(i+j,:), vol(:,:,i+j ), 'linear', 0));
            added = added + 1;
        end;
    end
end
end
