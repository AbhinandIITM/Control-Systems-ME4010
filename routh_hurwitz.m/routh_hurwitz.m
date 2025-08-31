%% Routh-Hurwitz Stability Criterion Table Generator V1.1

% Author: Varun Bhatia [bhatiav16@berkeley.edu]

% Date: July 28, 2019

% This function goes through the process of setting up a Routh-Hurwitz
% table to determine information regarding the in/stability of a control
% system given a closed/open-loop transfer function. More specifically, this
% will solve for and output how many closed/open-loop poles are in the
% right-half plane, the left-half plane, and on the jw-axis. Additionally,
% this function will take into account and generate the Routh_Hurwitz table
% given two special cases:
    % 1)First element of a row is 0. 
    % 2)Entire Row is 0.

%Recall that this function operates under the assumption that you are 
%inputting a CLOSED LOOP TRANSFER FUNCTION. Later developments will
%incorporate inputs of open loop forward transfer functions with unity
%and/or nonunity feedback.

%Any usage of MATLAB functions is accredited to MathWorks.

%UPDATE 1.1: This program now takes into account UNITY feedback.

%%

function []= routh_hurwitz()

prompt_n = 'Enter highest power in numerator: ';
num_hp = input(prompt_n);
num_hp_added = num_hp + 1;
num = zeros(1, num_hp_added);

for i = 1:num_hp+1
    num(i) = input('Enter coefficients, starting with the highest power: ');
end

prompt_d = 'Enter highest power in denominator: ';
den_hp = input(prompt_d);
den_hp_added = den_hp + 1;
den = zeros(1, den_hp_added);

for i = 1:den_hp+1
    den(i) = input('Enter coefficients, starting with the highest power: ');
end

prompt_TF = 'Are you entering an open loop or closed loop TF? If closed loop type "C", if open loop type "O" with quotes around letter: ';
TF_type = input(prompt_TF);

if TF_type == "O"
    
    H = tf(num,den);
    K = 1;
    G = feedback(H,K)
    [~,d] = tfdata(G, 'v');
    den = d; %to acquire the new denominator coefficients from post-feedback
    
elseif TF_type == "C"
    
    G = tf(num,den)
    
end



if den_hp == 0
    cols = 1;
elseif den_hp == 1 | den_hp == 2
    cols = 2;
elseif den_hp == 3 | den_hp == 4
    cols = 3;
elseif den_hp == 5 | den_hp == 6
    cols = 4;
elseif den_hp == 7 | den_hp == 8
    cols = 5;    
elseif den_hp == 9 | den_hp == 10
    cols = 6;
elseif den_hp == 11 | den_hp == 12
    cols = 7;
elseif den_hp == 13 | den_hp == 14
    cols = 8;  
end

%the matrix with proper number of rows and columns
RH_Table = zeros(den_hp+1, cols);


%set up the first element in row 1 & row 2 since MATLAB doesn't do 0 indexing
if cols == 1
    RH_Table(1,1) = den(1);
    
elseif cols > 1
    RH_Table(1,1) = den(1);
    RH_Table(2,1) = den(2);
    
end

if cols == 1
    fprintf('\n *Only one element in the RH Table, thus, you have a stable system.* \n');
    
elseif cols > 1
    first_count = 3; %since added first element (element 1) above so now skip one and start at 3rd element for first row
    second_count = 4; %since added second element (element 2) above so now skip one and start at 4th element for second row
    
     %filling in the first row of the RH table
    for i = 2:cols %start at 2 because filled in first element above
        if first_count <= den_hp_added
            RH_Table(1,i) = den(first_count);
            first_count = first_count + 2; 
        else 
            RH_Table(1,i) = 0;
        end
    end

    %filling in the second row of the RH table
    for i = 2:cols %start at 2 because filled in first element above
        if second_count <= den_hp_added
            RH_Table(2,i) = den(second_count);
            second_count = second_count + 2; 
        else 
            RH_Table(2,i) = 0;
        end
    end
    
    %now fill in remaining elements of RH Table with determinants
    X = zeros(2,2);
    
    for i = 1:(den_hp-1)
        for j = 1:(cols-1)
            X(1,1) = RH_Table(i,1);
            X(2,1) = RH_Table(i+1,1);
            X(1,2) = RH_Table(i,j+1);
            X(2,2) = RH_Table(i+1,j+1);
            
            RH_Table(i+2,j) = -det(X)/RH_Table(i+1,1);     
        end 
    end
    
    RH_Table % to show full RH table
    
    %now analyze the first column of the table for sign changes to
    %determine how many poles in the right half-plane.
    sign_changes = 0;
    for i = 1:(den_hp) 
        if (RH_Table(i,1) * RH_Table(i+1,1) < 0)
            sign_changes = sign_changes + 1;            
        end        
    end
        
    if sign_changes > 0 
         fprintf(['\n *You have %d right-half poles and %d left-half poles.',...
         ' Thus, you have an unstable system.* \n'],sign_changes,den_hp-sign_changes)
    else
        fprintf(['\n *You have %d right-half poles and %d left-half poles.',...
         ' Thus, you have a stable system.* \n'],sign_changes,den_hp-sign_changes) 
    end
        
end


%Special Case #1: Check if there exists a full row of zeros. If so, the
%following code will be utilized.

%first check if any array contains a row of zeros
row_of_zeros = 0;
tol = 1.e-6;
for i = 1:(den_hp + 1)
    if abs(RH_Table(i,:) - 0) < tol 
        row_of_zeros = 1;
        row = i;
        break; %break because you have identified a row of zeros
    end   
end


%If the check shows that you do indeed have a row of zeros, the following
%code will run

if (row_of_zeros == 1)
    row_above = RH_Table(row-1,:);
    aux_poly = poly2sym(row_above);
    diff_poly = diff(aux_poly);
    diff_poly_coeffs = flip(coeffs(diff_poly));
    
    for j = 1:cols
        if numel(diff_poly_coeffs) < cols
            diff_poly_coeffs = [diff_poly_coeffs, 0]; %add zero to the end of coefficients to make sure dimensions match when adding back to table
        else
            break;
        end
    end
    
    RH_Table(row,:) = diff_poly_coeffs; %replacing row in RH Table with the coefficients of row above, differentiated
     
    X = zeros(2,2);
    
    for i = row-1:(den_hp-1)
        for j = 1:(cols-1)
            X(1,1) = RH_Table(i,1);
            X(2,1) = RH_Table(i+1,1);
            X(1,2) = RH_Table(i,j+1);
            X(2,2) = RH_Table(i+1,j+1);
            
            RH_Table(i+2,j) = -det(X)/RH_Table(i+1,1);  
        end 
    end
    
    %Split up the sign changes process into two:
    %1) Check sign changes until two rows above row of zeros normally and
    %store these sign changes.
    %2) Check sign changes from one row above row of zeros down to bottom.
    %By symmetry, the right half poles must equal the left half poles, and
    %anything left over are jw-axis poles. Additionally, the row above the
    %row of zeros is required to be an even polynomial.
    
    %1)
    first_sign_changes = 0;
    second_sign_changes = 0;
    for k = 1:(row-2)
        if (RH_Table(k,1) * RH_Table(k+1,1) < 0)
            first_sign_changes = first_sign_changes + 1;
        end
    end
    
    %2)
    for k = row-1:den_hp
        if (RH_Table(k,1) * RH_Table(k+1,1) < 0)
            second_sign_changes = second_sign_changes + 1;
        end
    end
    
    total_sign_changes = first_sign_changes + second_sign_changes;
    left_half_poles  =((row-2)-first_sign_changes + second_sign_changes);
    jw_axis_poles = den_hp - (row-2) - (2*second_sign_changes);

    RH_Table
    
    if (total_sign_changes)  > 0
         fprintf(['\n *Upon changes made to the RH Table, you now have %d right-half',... 
             ' poles and %d left-half poles. Additionally, you have %d jw-axis poles. Thus, you have an unstable',...
             ' system.* \n'],total_sign_changes,left_half_poles,jw_axis_poles);
    elseif (total_sign_changes == 0)
        fprintf(['\n *Upon changes made to the RH Table, you now have %d',...
            ' right-half poles and %d left-half poles. Additionally, you have %d jw-axis poles.',...
         ' Thus, you have a marginally stable system.* \n'],total_sign_changes,left_half_poles, jw_axis_poles);
    end 
end


%Special Case #2: Check if zero exists in first column and if so, form reciprocal
%polynomial and redo process above to get determinants. Have a separate
%method for this portion, however, for the purposes of File Exchange, I've
%attached that portion below.

for i = 1:(den_hp+1)
    if RH_Table(i,1) == 0
        den = flip(den);
        if cols == 1
            RH_Table(1,1) = den(1);
        elseif cols > 1
            RH_Table(1,1) = den(1);
            RH_Table(2,1) = den(2);
        end
        
        if cols == 1
            fprintf('\n *Only one element in the RH Table, thus, you have a stable system.* \n');
        elseif cols > 1
            first_count = 3; 
            second_count = 4;
            for k = 2:cols 
                if first_count <= den_hp_added
                    RH_Table(1,k) = den(first_count);
                    first_count = first_count + 2; 
                else 
                    RH_Table(1,k) = 0;
                end
            end
            
            for l = 2:cols
                if second_count <= den_hp_added
                    RH_Table(2,l) = den(second_count);
                    second_count = second_count + 2; 
                else 
                    RH_Table(2,l) = 0;
                end
            end
            
            X = zeros(2,2);
            for m = 1:(den_hp-1)
                for j = 1:(cols-1)
                    X(1,1) = RH_Table(m,1);
                    X(2,1) = RH_Table(m+1,1);
                    X(1,2) = RH_Table(m,j+1);
                    X(2,2) = RH_Table(m+1,j+1);

                    RH_Table(m+2,j) = -det(X)/RH_Table(m+1,1);
                end 
            end
            
            sign_changes = 0;
            for k = 1:(den_hp)
                if (RH_Table(k,1) * RH_Table(k+1,1) < 0) %if the element * next element is negative, must be a sign change
                    sign_changes = sign_changes + 1;
                end
            end

            RH_Table
            if sign_changes > 0
                 fprintf(['\n *Upon changes made to the RH Table, you now have %d right-half',... 
                     ' poles and %d left-half poles. Thus, you have an unstable',...
                     ' system.* \n'],sign_changes,den_hp-sign_changes)
            else
                fprintf(['\n *Upon changes made to the RH Table, you now have %d',...
                    ' right-half poles and %d left-half poles.',...
                 ' Thus, you have a stable system.* \n'],sign_changes,den_hp-sign_changes)
            end 
        end
    end
end %end special case #2%


end


















