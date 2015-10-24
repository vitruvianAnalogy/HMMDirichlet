function [like_ratio_mean, like_ratio_12, like_ratio_21 ] = Comparison( mean_prob1,mean_prob2)
%HMM
[init_prob1,init_prob2,tran_prob1,tran_prob2,a_estimate_1,a_estimate_2,b_estimate_1,b_estimate_2,inferred_mean_pi_1,inferred_mean_pi_2,Final_inferred_mean_pi_1,Final_inferred_mean_pi_2,final_a1_inf,final_a2_inf,S1,S2]=forward_backward(mean_prob1,mean_prob2)

init_prob1;
tran_prob1;
%like_ratio_mean(30,10,2)=0;
like_ratio_mean=ones(30,10,2);
%like_ratio_12(30,10)=0;
like_ratio_12=ones(30,10);
%like_ratio_21(30,10)=0;
like_ratio_21=ones(30,10);
%For num of Zikj
%prob_obs1(30,10)=0;
%prob_obs2(30,10)=0;
prob_obs1=ones(30,10);
prob_obs2=ones(30,10);
%For mean HMM
%prob_obs_mean1(30,10)=0;
%prob_obs_mean2(30,10)=0;
prob_obs_mean1=ones(30,10);
prob_obs_mean2=ones(30,10);
%For R12
prob_obs12=ones(30,10);
%prob_obs12(30,10)=0;
%For R21
prob_obs21=ones(30,10);
%prob_obs21(30,10)=0;

for i=1:10
    for j=1:30
        prob_obs1(j,i)=1;
        prob_obs2(j,i)=1;
        prob_obs12(j,i)=1;
        prob_obs21(j,i)=1;
        %prob_obs121(j,i)=1;
        %prob_obs221(j,i)=1;
        %prob_obs_mean1(j,i)=1;
        %prob_obs_mean2(j,i)=1;
    end
end

for j=1:10    
        for l=1:30        
            for k=1:100
            %individual observation sequence representation
            %S1(l,:,j)
            %individual HMM
            %init_prob1(1,j)   P(S1) for jth HMM
            %tran_prob1(:,:,j)    P()
        
            if k~=1
                prob_obs1(l,j)= prob_obs1(l,j)*tran_prob1(S1(l,k,j),S1(l,k-1,j));
                prob_obs12(l,j)=prob_obs12(l,j)*final_a2_inf(S1(l,k,j),S1(l,k-1,j));
                prob_obs21(l,j)=prob_obs21(l,j)*final_a1_inf(S2(l,k,j),S2(l,k-1,j));
                prob_obs2(l,j)= prob_obs2(l,j)*tran_prob2(S2(l,k,j),S2(l,k-1,j));
                prob_obs_mean1(l,j)=prob_obs_mean1(l,j)*final_a1_inf(S1(l,k,j),S1(l,k-1,j));
                prob_obs_mean2(l,j)=prob_obs_mean2(l,j)*final_a2_inf(S2(l,k,j),S2(l,k-1,j));
            end
            prob_obs1(l,j)=prob_obs1(l,j)*init_prob1(S1(l,1,j),j);
            prob_obs2(l,j)=prob_obs2(l,j)*init_prob2(S2(l,1,j),j);            
            %prob_obs12(l,j)=prob_obs12(l,j)*Final_inferred_mean_pi_2(S1(l,1,j),j);
            %prob_obs21(l,j)=prob_obs21(l,j)*Final_inferred_mean_pi_1(S2(l,1,j),j);
            prob_obs12(l,j)=prob_obs12(l,j)*Final_inferred_mean_pi_2(S1(l,1,j));
            prob_obs21(l,j)=prob_obs21(l,j)*Final_inferred_mean_pi_1(S2(l,1,j));
            %prob_obs_mean1(l,j)=prob_obs_mean1(l,j)*Final_inferred_mean_pi_1(S1(l,1,j),j);
            %prob_obs_mean2(l,j)=prob_obs_mean2(l,j)*Final_inferred_mean_pi_2(S2(l,1,j),j);
            prob_obs_mean1(l,j)=prob_obs_mean1(l,j)*Final_inferred_mean_pi_1(S1(l,1,j));
            prob_obs_mean2(l,j)=prob_obs_mean2(l,j)*Final_inferred_mean_pi_2(S2(l,1,j));
            end
            
        end
end
for j=1:10
    for l=1:30
        like_ratio_mean(l,j,1)=prob_obs1(l,j)/prob_obs_mean1(l,j)
        like_ratio_mean(l,j,2)=prob_obs2(l,j)/prob_obs_mean2(l,j)
        like_ratio_12(l,j)=prob_obs_mean1(l,j)/prob_obs12(l,j)
        like_ratio_21(l,j)=prob_obs_mean2(l,j)/prob_obs21(l,j)
    end
end

k=0;
figure(1);
for i=1:10
    for j=1:2
        k=k+1;
        subplot(4,5,k)
        semilogy(like_ratio_mean(:,i,j));        
    end    
end
figure(2);
k=0;
for i=1:10
    k=k+1;
    subplot(3,4,k)    
    semilogy(like_ratio_12(:,i));    
end
figure(3);
k=0;
for i=1:10
    k=k+1;
    subplot(2,5,k)    
    semilogy(like_ratio_21(:,i));
end


