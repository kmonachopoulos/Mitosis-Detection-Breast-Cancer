%-------------------------------------------------------------------------% 
%  Project       : Mitosis-Detection-Breast-Cancer                        %
%  File          : threshold_of iterations_vs_otsu.m                      %
%  Description   : threshold by iterations                                %
%  Author        : Monachopoulos Konstantinos                              %
%-------------------------------------------------------------------------%

function [ Output_Image ] = thresh_by_iterations( Input_Image )

%% Initialization Variables

Init_Image=Input_Image(:,:,1);
Histogram_Init_Image=myhist_fcn(Init_Image);
Cumulative_Histogram_Init_Image=mysumhist_fcn(Histogram_Init_Image);
k=1:256;

% Choose of initial value of threshold (etc. mean value)

% We multiply every pixel with the luminance allocated to them /  whole of pixels
Threshold_by_iterations=round(sum(Histogram_Init_Image.*k)/Cumulative_Histogram_Init_Image(end));

for i=1:2
    
    
    %% classification of pixels into two groups R1 and R2 depending of T threshold  
	
    R1_Histogram=Histogram_Init_Image(1:Threshold_by_iterations);
    R1_Cumulative=mysumhist_fcn(R1_Histogram);
    R2_Histogram=Histogram_Init_Image(Threshold_by_iterations+1:length(Histogram_Init_Image));
    R2_Cumulative=mysumhist_fcn(R2_Histogram);
    
    %% Calculation of mean values ì1 and ì2 of the two groups
    
    m1(i)=sum(R1_Histogram.*[1:Threshold_by_iterations])/R1_Cumulative(end);
    m2(i)=sum(R2_Histogram.*[Threshold_by_iterations+1:256])/R2_Cumulative(end);
    
    %% selection of new threshold value Ô= 1/2 (ì1+ì2)
    
    Threshold_by_iterations=round((m1(i)+m2(i))/2);
    
end

while(m1(i)~=m1(i-1) && m2(i)~= m2(i-1))
    
    i=i+1;
    R1_Histogram=Histogram_Init_Image(1:Threshold_by_iterations);
    R1_Cumulative=mysumhist_fcn(R1_Histogram);
    R2_Histogram=Histogram_Init_Image(Threshold_by_iterations+1:length(Histogram_Init_Image));
    R2_Cumulative=mysumhist_fcn(R2_Histogram);
    
    m1(i)=sum(R1_Histogram.*[1:Threshold_by_iterations])/R1_Cumulative(end);
    m2(i)=sum(R2_Histogram.*[Threshold_by_iterations+1:256])/R2_Cumulative(end);
    
    Threshold_by_iterations=round((m1(i)+m2(i))/2);
end

Threshold_by_iterations_Image=Init_Image;
Threshold_by_iterations_Image(Threshold_by_iterations_Image>=Threshold_by_iterations)=255;
Threshold_by_iterations_Image(Threshold_by_iterations_Image<Threshold_by_iterations)=0;
Output_Image=~Threshold_by_iterations_Image;