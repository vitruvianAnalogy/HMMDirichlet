function [init_prob1,init_prob2,tran_prob1,tran_prob2,a_estimate_1,a_estimate_2,b_estimate_1,b_estimate_2,inferred_mean_pi_1,inferred_mean_pi_2,Final_inferred_mean_pi_1,Final_inferred_mean_pi_2,final_a1_inf,final_a2_inf,S1,S2] = forward_backward( mean_prob1,mean_prob2)
%forward_backward:estimates transition probabilities and initial
%probabilities
S1=zeros(30,100,10);
%observation sequence mean 1
O1=zeros(30,100,10);
S2=zeros(30,100,10);
%observation sequence mean 2
O2=zeros(30,100,10);
%Generate 20 HMM's
for i=1:10
    [S1(:,:,i),O1(:,:,i),init_prob1(:,i),tran_prob1(:,:,i)]=generate_state_obs(mean_prob1);
    [S2(:,:,i),O2(:,:,i),init_prob2(:,i),tran_prob2(:,:,i)]=generate_state_obs(mean_prob2);
end

len1=size(S1);
%len2=size(S2);
eps1=zeros(99,5,5,10);
eps2=zeros(99,5,5,10);
%Initialise epsilon
%First define eps1(t,i,j) as the probability of being in state i at time t and state i+1 in time t+1
for num_hmm=1:10
  for t=1:(len1(2)-1)      %upto 99
      for num_seq=1:len1(1) %upto 30
         for i=1:5    %upto 5
           for j=1:5  %upto 5
            if (S1(num_seq,t,num_hmm)==i)&&(S1(num_seq,t+1,num_hmm)==j)
                eps1(t,i,j,num_hmm)=eps1(t,i,j,num_hmm)+1;
            end
            if (S2(num_seq,t,num_hmm)==i)&&(S2(num_seq,t+1,num_hmm)==j)
                eps2(t,i,j,num_hmm)=eps2(t,i,j,num_hmm)+1;
            end
           end
         end
      end
  end
end

%Initialise gamma array
Gamma1_a=zeros(99,5);
Gamma2_a=zeros(99,5);
for num_hmm=1:10  
    for t=1:(len1(2)-1)   %upto 100
       for i=1:5
         for j=1:5
             Gamma1_a(t,i)=Gamma1_a(t,i)+eps1(t,i,j,num_hmm);
             Gamma2_a(t,i)=Gamma2_a(t,i)+eps2(t,i,j,num_hmm);
         end
       end
    end
end
%%%Init_mean part-Qn4
Init_Gamma1_a=zeros(99,5,10);
Init_Gamma2_a=zeros(99,5,10);

for num_hmm=1:10
     for t=1:(len1(2)-1)   %upto 100
       for i=1:5
         for j=1:5
             Init_Gamma1_a(t,i,num_hmm)=Init_Gamma1_a(t,i,num_hmm)+eps1(t,i,j,num_hmm);
             Init_Gamma2_a(t,i,num_hmm)=Init_Gamma2_a(t,i,num_hmm)+eps2(t,i,j,num_hmm);
         end
       end
    end
end

for num_hmm=1:10
    for i=1:5
        inferred_mean_pi_1(i,num_hmm)= Init_Gamma1_a(1,i,num_hmm);
        inferred_mean_pi_2(i,num_hmm)= Init_Gamma2_a(1,i,num_hmm);
    end
end
Final_inferred_mean_pi_1=mean(inferred_mean_pi_1,2)/30;
Final_inferred_mean_pi_2=mean(inferred_mean_pi_2,2)/30;
inferred_mean_pi_1=inferred_mean_pi_1/30;
inferred_mean_pi_2=inferred_mean_pi_2/30;
for num_hmm=1:10
  for i=1:5
       sum_gamma1_inf(num_hmm)=0;
       sum_gamma2_inf(num_hmm)=0;
       for t=1:99
        sum_gamma1_inf(num_hmm)=sum_gamma1_inf(num_hmm)+Init_Gamma1_a(t,i,num_hmm);
        sum_gamma2_inf(num_hmm)=sum_gamma2_inf(num_hmm)+Init_Gamma2_a(t,i,num_hmm);    
        end
         for j=1:5
          sum_eps1_inf(num_hmm)=0;
          sum_eps2_inf(num_hmm)=0;
            for t=1:99
            sum_eps1_inf(num_hmm)=sum_eps1_inf(num_hmm)+eps1(t,i,j,num_hmm);
            sum_eps2_inf(num_hmm)=sum_eps2_inf(num_hmm)+eps2(t,i,j,num_hmm);
            end
            a1_inf(i,j,num_hmm)=sum_eps1_inf(num_hmm)/sum_gamma1_inf(num_hmm);
            a2_inf(i,j,num_hmm)=sum_eps2_inf(num_hmm)/sum_gamma2_inf(num_hmm);
          end
    end
end
final_a1_inf=mean(a1_inf,3);
final_a2_inf=mean(a1_inf,3);

%Initialize gamma array for b
Gamma1_b=zeros(5,100,10);
Gamma2_b=zeros(5,100,10);
for num_hmm=1:10
  for num_seq=1:30
    for t=1:100
        temp_state1=S1(num_seq,t,num_hmm);
        temp_obs1=O1(num_seq,t,num_hmm);
        temp_state2=S2(num_seq,t,num_hmm);
        temp_obs2=O2(num_seq,t,num_hmm);
        for curr_indx=1:100
           if (S1(num_seq,curr_indx,num_hmm)==temp_state1) && (O1(num_seq,curr_indx,num_hmm)==temp_obs1)
                       Gamma1_b(temp_state1,t,num_hmm)=Gamma1_b(temp_state1,t,num_hmm)+1;
           end
           if (S2(num_seq,curr_indx,num_hmm)==temp_state2) && (O2(num_seq,curr_indx,num_hmm)==temp_obs2)
                       Gamma2_b(temp_state2,t,num_hmm)=Gamma2_b(temp_state2,t,num_hmm)+1;
           end
        end
     end
   end
end
               
%Compute A
a1=zeros(5,5,10);
a2=zeros(5,5,10);
b1=zeros(5,10);
b2=zeros(5,10);

for num_hmm=1:10
  for i=1:5
       sum_gamma1=0;
       sum_gamma2=0;
       sum_gamma1_b=0;
       sum_gamma2_b=0;
        for t=1:99
        sum_gamma1=sum_gamma1+Gamma1_a(t,i);
        sum_gamma2=sum_gamma2+Gamma2_a(t,i);    
        end
         for t=1:100
             sum_gamma1_b=sum_gamma1_b+Gamma1_b(i,t,num_hmm);
             sum_gamma2_b=sum_gamma2_b+Gamma2_b(i,t,num_hmm);
         end
          for j=1:5
          sum_eps1=0;
          sum_eps2=0;
            for t=1:99
            sum_eps1=sum_eps1+eps1(t,i,j,num_hmm);
            sum_eps2=sum_eps2+eps2(t,i,j,num_hmm);
            end
            a1(i,j,num_hmm)=sum_eps1/sum_gamma1;
            a2(i,j,num_hmm)=sum_eps2/sum_gamma2;
          end
          b1(i,num_hmm)=sum_gamma1_b/sum_gamma1;
          b2(i,num_hmm)=sum_gamma2_b/sum_gamma2;
  end
end
init_a1=a1;
init_a2=a2;
init_b1=b1;
init_b2=b2;
temp_a1=zeros(5,5,10);
temp_a2=zeros(5,5,10);
temp_b1=zeros(5,10);
temp_b2=zeros(5,10);

max_diff_b1=0.5;
max_diff_b2=0.5;
max_diff_a1=0.5;
max_diff_a2=0.5;
iteration_count=0;

while max_diff_b1>0.05 || max_diff_b2>0.05 || max_diff_a1>0.05 || max_diff_a2>0.05 
    temp_a1=a1;
    temp_a2=a2;
    temp_b1=b1;
    temp_b2=b2;
%Find alpha
alpha_S1=zeros(5,100,10);
alpha_S2=zeros(5,100,10);

for num_hmm=1:10
  for i=1:5
    alpha_S1(i,1,num_hmm)=init_prob1(i,num_hmm)*b1(i,num_hmm);
    alpha_S2(i,1,num_hmm)=init_prob2(i,num_hmm)*b2(i,num_hmm);
  end
end
for num_hmm=1:10
  for t=1:99
    for j=1:5
        prod_sum1=0;
        prod_sum2=0;
        for i=1:5
             prod_sum1=prod_sum1+(alpha_S1(i,t,num_hmm)*a1(i,j,num_hmm));
             prod_sum2=prod_sum2+(alpha_S2(i,t,num_hmm)*a2(i,j,num_hmm));
        end
        alpha_S1(j,t+1,num_hmm)=prod_sum1*b1(j,num_hmm);
        alpha_S2(j,t+1,num_hmm)=prod_sum2*b2(j,num_hmm);
    end
  end
end
prob_O_lambda_1=zeros(10);
prob_O_lambda_2=zeros(10);
for num_hmm=1:10
    for i=1:5 
      prob_O_lambda_1(num_hmm)=prob_O_lambda_1(num_hmm)+alpha_S1(i,100,num_hmm);
      prob_O_lambda_2(num_hmm)=prob_O_lambda_2(num_hmm)+alpha_S2(i,100,num_hmm);
    end
end
prob_O_lambda_1=prob_O_lambda_1(1:10,1:1);
prob_O_lambda_2=prob_O_lambda_2(1:10,1:1);
%Find beta
beta1=zeros(5,100);
beta2=zeros(5,100);
sum_a1_b1=0;
sum_a2_b2=0;
for num_hmm=1:10
  for i=1:5
    beta1(i,100,num_hmm)=1;
    beta2(i,100,num_hmm)=1;
  end
end
for num_hmm=1:10
 for t=99:-1:1
    for i=1:5
     for j=1:5
        beta1(i,t,num_hmm)=beta1(i,t,num_hmm)+(a1(i,j,num_hmm)*b1(j,num_hmm));
        beta2(i,t,num_hmm)=beta2(i,t,num_hmm)+(a2(i,j,num_hmm)*b2(j,num_hmm));
     end
    end
 end
end
for num_hmm=1:10
  for i=1:5
    for j=1:5
        for t=2:100
            Gamma_r_1(i,j,t,num_hmm)=(alpha_S1(i,t-1,num_hmm)*a1(i,j,num_hmm)*b1(j,num_hmm)*beta1(j,t,num_hmm))/prob_O_lambda_1(num_hmm);
            Gamma_r_2(i,j,t,num_hmm)=(alpha_S2(i,t-1,num_hmm)*a2(i,j,num_hmm)*b2(j,num_hmm)*beta2(j,t,num_hmm))/prob_O_lambda_2(num_hmm);
        end
    end
  end
end

%Initialize gamma array for b
%Initialize gamma array for b
Gamma1_b=zeros(5,100,10);
Gamma2_b=zeros(5,100,10);
for num_hmm=1:10
  for num_seq=1:30
    for t=1:100
        temp_state1=S1(num_seq,t,num_hmm);
        temp_obs1=O1(num_seq,t,num_hmm);
        temp_state2=S2(num_seq,t,num_hmm);
        temp_obs2=O2(num_seq,t,num_hmm);
        for curr_indx=1:100
           if (S1(num_seq,curr_indx,num_hmm)==temp_state1) && (O1(num_seq,curr_indx,num_hmm)==temp_obs1)
               for j=1:5
                       Gamma1_b(temp_state1,t,num_hmm)=Gamma_r_1(temp_state1,j,t,num_hmm)+1;
               end
           end
           if (S2(num_seq,curr_indx,num_hmm)==temp_state2) && (O2(num_seq,curr_indx,num_hmm)==temp_obs2)
               for j=1:5        
                       Gamma2_b(temp_state2,t,num_hmm)=Gamma_r_2(temp_state2,j,t,num_hmm)+1;
               end
           end
        end
     end
   end
end
for num_hmm=1:10
    for i=1:5
        sum_gamma1_b=0;
        sum_gamma2_b=0;
      
         for t=1:100
             sum_gamma1_b=sum_gamma1_b+Gamma1_b(i,t,num_hmm);
             sum_gamma2_b=sum_gamma2_b+Gamma2_b(i,t,num_hmm);
         end
        for j=1:5
        sum_gamma_1_r=0;
        sum_gamma_2_r=0;
        sum_gamma_1_dnr=0;
        sum_gamma_2_dnr=0;
            for t=2:100
                sum_gamma_1_r=sum_gamma_1_r+Gamma_r_1(i,j,t,num_hmm);
                sum_gamma_2_r=sum_gamma_2_r+Gamma_r_2(i,j,t,num_hmm);
                
                for k=1:5
                sum_gamma_1_dnr=sum_gamma_1_dnr+Gamma_r_1(i,k,t,num_hmm);
                sum_gamma_2_dnr=sum_gamma_2_dnr+Gamma_r_2(i,k,t,num_hmm);   
                end
            end
            
            a_estimate_1(i,j,num_hmm)=sum_gamma_1_r/sum_gamma_1_dnr;
            a_estimate_2(i,j,num_hmm)=sum_gamma_2_r/sum_gamma_2_dnr;
            end
            b_estimate_1(i,num_hmm)=sum_gamma1_b/sum_gamma_1_dnr;
            b_estimate_2(i,num_hmm)=sum_gamma2_b/sum_gamma_2_dnr;
        
    end
end
a1=a_estimate_1;
a2=a_estimate_2;
b1=b_estimate_1;
b2=b_estimate_2;
diff_a1=(a1-temp_a1);
diff_a2=(a2-temp_a2);
diff_b1=(b1-temp_b1);
diff_b2=(b2-temp_b2);
max_diff_a1_vector=max(diff_a1);
max_diff_a2_vector=max(diff_a2);
max_diff_b1_vector=max(diff_b1);
max_diff_b2_vector=max(diff_b2);
max_diff_a1_array=max(max_diff_a1_vector);
max_diff_a2_array=max(max_diff_a2_vector);
max_diff_b1=max(max_diff_b1_vector)
max_diff_b2=max(max_diff_b2_vector)
max_diff_a1=max(max_diff_a1_array)
max_diff_a2=max(max_diff_a2_array)

max_diff_b1=abs(max_diff_b1)
max_diff_b2=abs(max_diff_b2)
max_diff_a1=abs(max_diff_a1)
max_diff_a2=abs(max_diff_a2)

            
 iteration_count=iteration_count+1           
 
            
            
end