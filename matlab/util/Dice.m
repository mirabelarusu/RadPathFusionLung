function D = Dice( argins )
%Dice Returns the Dice coefficient of the data 
%   Dice(A,B) = 2*JaccardIndex(A,B)/(1+JaccardIndex(A,B))
%   
% MR


dims  = size(size(argins),2);

if (dims<4)
    disp(['Input must be either a 4D dataset or two 3D datasets ']);
    JI = 0;
    return;
end

data = argins;

I = data(:,:,:,1);
s = sumall(I);

for i = 2:size(argins,dims)
    I = I & data(:,:,:,i);
    s = s + sumall(data(:,:,:,i));
end

if s ~= 0
    D = size(data,4) *sumall(I)/s;
else
    D = 0;
end

end

