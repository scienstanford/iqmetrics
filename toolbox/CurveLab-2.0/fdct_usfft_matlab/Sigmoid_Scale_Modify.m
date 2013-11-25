function En_Co = Sigmoid_Scale_Modify(Coeff,N_Coeff,p,s,c,type)

%% 单层尺度系数调整
switch type
    case {'Mu_L'}
        En_Co = cell(1, length(Coeff));
        En_Co{1} = Coeff{1};
        for i = 2:length(Coeff)
            for j = 1:length(Coeff{i})
                En_Co{i}{j} = (Coeff{i}{j}./abs(Coeff{i}{j})).*abs(Coeff{i}{j}).^p;
            end
        end 
    case {'Mu_S'}
        En_Co = cell(1, length(Coeff));
        En_Co{1} = Coeff{1};        
        for i = 2:length(Coeff)
            maxp = max(Coeff{i}(:));            
            En_Co{i} = Coeff{i};
            [xray,yray] = find(abs(Coeff{i})<=maxp);
            index_num = length(xray);
            while index_num~=0
                if i<4                
                    En_Co{i}(xray(index_num),yray(index_num)) = 3*Coeff{i}(xray(index_num),yray(index_num))*((1-abs(Coeff{i}(xray(index_num),yray(index_num)))/maxp)^p)+Coeff{i}(xray(index_num),yray(index_num));
                else
                    En_Co{i} = 3*Coeff{i}(xray(index_num),yray(index_num))*((1-abs(Coeff{i}(xray(index_num),yray(index_num)))/maxp)^p)+Coeff{i}(xray(index_num),yray(index_num));
                end
                index_num = index_num-1;
            end
        end
    case {'Mu_Cur'}        
        En_Co = cell(1, length(Coeff));
        En_Co{1} = Coeff{1};
        for i = 1:length(Coeff)
            for j = 1:length(Coeff{i})
                maxp = max(abs(Coeff{i}{j}(:)));
                comc = maxp*0.5; %30; 
                C_delta = c*N_Coeff{i}{j};
                En_Co{i}{j} = Coeff{i}{j};               
                [xray,yray] = find(abs(Coeff{i}{j})>=C_delta&abs(Coeff{i}{j})<2*C_delta);
                index_num = length(xray);
                while index_num~=0
                    Yc = (abs(Coeff{i}{j}(xray(index_num),yray(index_num)))-C_delta)/C_delta*(comc/C_delta)^p+(2*C_delta-abs(Coeff{i}{j}(xray(index_num),yray(index_num))))/C_delta;
                    En_Co{i}{j}(xray(index_num),yray(index_num)) = Coeff{i}{j}(xray(index_num),yray(index_num))*Yc;
                    index_num = index_num-1;
                end
                [xray,yray] = find(abs(Coeff{i}{j})<comc&abs(Coeff{i}{j})>=2*C_delta);
                index_num = length(xray);
                while index_num~=0
                    Yc = (comc/abs(Coeff{i}{j}(xray(index_num),yray(index_num))))^p;
                    En_Co{i}{j}(xray(index_num),yray(index_num)) = Coeff{i}{j}(xray(index_num),yray(index_num))*Yc;
                    index_num = index_num-1;
                end  
                [xray,yray] = find(abs(Coeff{i}{j})>=comc);
                index_num = length(xray);
                while index_num~=0
                    Yc = (comc/abs(Coeff{i}{j}(xray(index_num),yray(index_num))))^s;
                    En_Co{i}{j}(xray(index_num),yray(index_num)) = Coeff{i}{j}(xray(index_num),yray(index_num))*Yc;
                    index_num = index_num-1;
                end  
            end
         end
          
end
