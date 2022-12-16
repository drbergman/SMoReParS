function I = buildIndices(ndims)


if ndims == 3
    I.subs = 1:3;
    I.ind = 4;
    I.event = 5;
    I.phase = 6;

else
    I.subs = 1:2;
    I.ind = 3;
    I.event = 4;
    I.phase = 5;
end
