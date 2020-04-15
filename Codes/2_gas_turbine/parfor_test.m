% 

% Test parfor

clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parfor tips
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Defining data to be used within parfor
% Matrices are not great. 
% Scalars and vectors are good. 

%%% Storing data within parfor
% Whatever you compute within parfor won't make it out unless you properly
% store it. 
% Use cell arrays - store data in a cell array that is as long as how many
% parfor loops there are. 
% If your parfor goes through n different loops, preallocate - don't follow
% that. 

%%% Orgainizing the data later
% If you need to look through all your data, like if each cell was some
% state ... simply loop through the length of myData and do what you need
% to. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Number to loop through with parfor

n = 10; 

%%% Data that I'll definitely need later; 

crucialData = 5; 
crucialVector = [1,2,3]; 

%%% Preallocating cell array "myData"
myData{n}=[]; 

%%% Looping through with parfor
parfor ii = 1:n
    myData{ii}.d = crucialData+ii; 
    
    myData{ii}.vec = [ii ii ii] + crucialVector
    
    for i_1 = 1:5
        i_1
    end
end

newData = zeros(n,1); 
for ii = 1:length(myData)
    newData(ii) = myData{ii}.d + sum(myData{ii}.vec); 
end
