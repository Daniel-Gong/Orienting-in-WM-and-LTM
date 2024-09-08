%Bootstrap Method compare two conditions
function[difference,p_val] = bootstrp_compare(data1,data2)
data1_resample = bootstrp(10000,@mean,data1);
data2_resample = bootstrp(10000,@mean,data2);
difference = data1_resample - data2_resample;
p_val = sum(difference > 0)/length(difference);
if p_val > 0.5
    p_val = 1 - p_val;
end
p_val = 2 * p_val;



    
