%-------------------------------------------------------------------------%
%  Project       : Mitosis-Detection-Breast-Cancer                        %
%  File          : mysumhist_fcn.m                                        %
%  Description   : Export of image cumulative histogram                   %
%  Author        : Monachopoulos Konstantinos                              %
%-------------------------------------------------------------------------%

function [ cumulative_hist ] = mysumhist_fcn( hist )
cumulative_hist(1)=hist(1);
for i=2:size(hist,2)
    cumulative_hist(i)=cumulative_hist(i-1)+hist(i);
end
end

