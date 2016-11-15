
function [d] = Averageslices(gradients, slice, d)
gradients = gradients+1;
slice = slice+1;

for x=1:size(d,1)
    for y=1:size(d,2)

        d(x,y,slice,gradients) = (d(x,y,slice-1,gradients) + d(x,y,slice+1,gradients))/2;
    end
end