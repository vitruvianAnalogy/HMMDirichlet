function [S,O,init_prob,tran_prob]= generate_state_obs(mean_prob)
%mean_prob=[.33 .11 .19 .21 .16];
    [tran_prob, init_prob] = draw_dirichlet_markov(mean_prob);
%    disp(init_prob);
%    disp(tran_prob);
    C=round(init_prob*30);
%    disp(C);
    
    for n=2:5
    C(n,1)=C(n-1,1)+C(n,1);
    end
    S1=zeros(30,1);
   S1(1:C(1))=1;

   for j=2:5
    S1(C(j-1)+1:C(j))=j;
   end
    S1=S1(1:30,1:1);
    S2=zeros(30,99);
    for j=1:30
        for k=1:99
            S2(j,k)=randi(5,1,1);
        end
    end
    S=horzcat(S1,S2);
    S(S==0)=randi(5,1,1);
    O=zeros(30,100);
    for j=1:30
        for k=1:100
            O(j,k)=normrnd(S(j,k),0.3);
        end
    end
%        disp(S1);
%        disp(S);
%        disp(O);
end
