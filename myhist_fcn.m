%-------------------------------------------------------------------------%
%  Project       : Mitosis-Detection-Breast-Cancer                        %
%  File          : my_hist_fcn.m                                          %
%  Description   : Export of image histogram                              %
%  Author        : Monahopoulos Konstantinos                              %
%-------------------------------------------------------------------------%

function [ hist ] = myhist_fcn(myimage )

for greylevel=0:255
    counter=0;
    for i=1:size(myimage,1)
        for j=1:size(myimage,2)
            
            if myimage(i,j)== greylevel;
                counter=counter+1;
            end
        end
    end
    hist(greylevel+1)=counter;
end

