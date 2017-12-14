function V = opti_flow_poly(I,BC,cond_fact)

[Y,X,T] = size(I);
D = [];
J = [];
V1 = zeros(Y,X,T);
V2 = zeros(Y,X,T);
 %% calculate D and J
for tt = -(round(T/2) - 1):(round(T/2) - 1)
    for yy = -(round(Y/2)-1):(round(Y/2) - 1)
        for xx = -(round(X/2)-1):(round(X/2)-1)
   
            D = [D; 1 , xx , yy , tt , (xx)^2 , xx*yy , (yy)^2 , yy*tt , (tt)^2 , xx*tt , (xx)^3 , ((xx)^2)*yy , xx*(yy)^2 , (yy)^3 , ((yy)^2)*tt , yy*(tt)^2 , (tt)^3 , ((xx)^2)*tt , xx*(tt)^2 , xx*yy*tt ];
            J = [ J; I( round(Y/2)-yy , round(X/2)+xx, round(T/2)+tt ) ];
            
        end
    end
end
%% get a1 - a20 coefficient values
AP = pinv(D)*double(J);

if( cond(AP) < cond_fact )%% if not set OF to zero
    
a1 = AP(1); a2 = AP(2); a3 = AP(3); a4 = AP(4); a5 = AP(5); a6 = AP(6); a7 = AP(7); a8 = AP(8); a9 = AP(9); a10 = AP(10); a11 = AP(11); a12 = AP(12); a13 = AP(13);a14 = AP(14); a15 = AP(15); a16 = AP(16); a17 = AP(17); a18 = AP(18); a19 = AP(19); a20 = AP(20);
%% get OF values
for tt = -(round(T/2) - 1):(round(T/2) - 1)
    for yy = -(round(Y/2)-1):(round(Y/2) - 1)
        for xx = -(round(X/2)-1):(round(X/2)-1)
            
        [MM,~] = find((sum(cell2mat(BC(:,end)) == [xx yy tt],2)) == 3);            

        A = eval(BC{MM,1});     %% eval values of placeholder variables with calculayted coefficients
        b = eval(BC{MM,2});     %% A and b calculated for the image gradients
              
        VV = pinv(A)*b;
        
        if(cond(A) < cond_fact )
        V1(round(Y/2)-yy , round(X/2)+xx, round(T/2)+tt) = VV(1);       %% if condition too high then OF zero
        V2(round(Y/2)-yy , round(X/2)+xx, round(T/2)+tt) = VV(2);
        end
        
        end
    end
end

end

V = cat(3,V1,V2);       %% V1 and V2

end