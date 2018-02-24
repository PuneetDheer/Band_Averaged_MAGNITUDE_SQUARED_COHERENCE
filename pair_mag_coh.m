% PAIRWISE COMPUTATION MAGNITUDE SQUARED COHERENCE (BAND_AVERAGED)
% Date:01-August-2017
%
% INPUT:
% channel_matrix= row wise data channel
% Ws = Window size (in sample point)
% Lp = shift the window by sample point
% FREQ_RANGE = frequency range of interest for average coherence like[30:80]
%
% OUTPUT:
% coherence_matrix = pair wise BAND_AVERAGED MAGNITUDE SQUARED COHERENCE
% windowed_eigen_vector = all eigen vector time series
% H_e_v = Highest eigen value for common synchronization
% S_esti_MI = S-estimator based synchronization 
%%

function [ coherence_matrix,windowed_eigen_vector, H_e_v, S_esti_coh ] = pair_mag_coh( channel_matrix,Ws,Lp,FREQ_RANGE)


n=size(channel_matrix,1);%no. of channels
windowed_eigen_vector=[];
coherence_matrix=[];
coherence_matrix_across=ones(n,n);
tic
numChannels = size(channel_matrix, 1);

for i = 1:numChannels
    
     comp(i,:) = channel_matrix(i, :);
    
end

for i=1:n
    for j=i+1:n
        
        [co]= BA_msc(comp(i,:),comp(j,:),Ws,Lp,FREQ_RANGE);
               
        coherence_matrix=[coherence_matrix;co];
           
    end
end

p=size(coherence_matrix,2); %no. of windows


for i=1:p %no. of windows
    a=0;
    for j=1:n
        for k=j+1:n
            a=a+1;
            coherence_matrix_across(j,k)=coherence_matrix(a,i);     
            coherence_matrix_across(k,j)=coherence_matrix(a,i);
            
        end
    end
    diag(coherence_matrix_across);%just to check
   %Note: here sum of eigen value should be equal to the no. of Channels.
   %sum of diagonal of windowed MI Matrix is also equal to the no. of Channels.
    windowed_matrix=coherence_matrix_across;
    coherence_matrix_across(isnan(coherence_matrix_across))=0;
%     pause
    V=sort(eig(coherence_matrix_across));
    windowed_eigen_vector=[windowed_eigen_vector V];
    coherence_matrix_across=ones(n,n);
    norm_eig = V./sum(V);
    S_esti_coh(i)=real(1+(sum(norm_eig.*log(norm_eig))/log(n)));
    
end

H_e_v=windowed_eigen_vector(end,:);

toc

end

