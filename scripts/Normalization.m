function se = Normalization(data)
    Num_Cond = size(data,2);
    Nsub = size(data,1);
    for nsub=1:Nsub
        sub_mean(nsub) = mean(data(nsub,:));
    end
    group_mean = mean(mean(data));
    for nsub = 1:Nsub
        for cond = 1:Num_Cond
            Norm(nsub,cond)=data(nsub,cond)-sub_mean(nsub)+group_mean;
        end
    end
    for cond = 1:Num_Cond
        se(cond) = std(Norm(:,cond))/sqrt(Nsub);
    end
end
% for nsub = 1:Nsub
%     sub_mean(nsub) = mean([MEAN_LOW_distance0(nsub) MEAN_LOW_distance1(nsub) MEAN_LOW_distance2(nsub) MEAN_LOW_distance3(nsub) MEAN_LOW_distance4(nsub) MEAN_MID_distance0(nsub) MEAN_MID_distance1(nsub) MEAN_MID_distance2(nsub) MEAN_MID_distance3(nsub) MEAN_MID_distance4(nsub) MEAN_HIGH_distance0(nsub) MEAN_HIGH_distance1(nsub) MEAN_HIGH_distance2(nsub) MEAN_HIGH_distance3(nsub) MEAN_HIGH_distance4(nsub)]);
% end
% group_mean = mean([mean(MEAN_LOW_distance0) mean(MEAN_LOW_distance1) mean(MEAN_LOW_distance2) mean(MEAN_LOW_distance3) mean(MEAN_LOW_distance4) mean(MEAN_MID_distance0) mean(MEAN_MID_distance1) mean(MEAN_MID_distance2) mean(MEAN_MID_distance3) mean(MEAN_MID_distance4) mean(MEAN_HIGH_distance0) mean(MEAN_HIGH_distance1) mean(MEAN_HIGH_distance2) mean(MEAN_HIGH_distance3) mean(MEAN_HIGH_distance4)]);
% for nsub = 1:Nsub
%     Norm_low0(nsub) = MEAN_LOW_distance0(nsub) - sub_mean(nsub) + group_mean;
%     Norm_low1(nsub) = MEAN_LOW_distance1(nsub) - sub_mean(nsub) + group_mean;
%     Norm_low2(nsub) = MEAN_LOW_distance2(nsub) - sub_mean(nsub) + group_mean;
%     Norm_low3(nsub) = MEAN_LOW_distance3(nsub) - sub_mean(nsub) + group_mean;
%     Norm_low4(nsub) = MEAN_LOW_distance4(nsub) - sub_mean(nsub) + group_mean;
%     Norm_mid0(nsub) = MEAN_MID_distance0(nsub) - sub_mean(nsub) + group_mean;
%     Norm_mid1(nsub) = MEAN_MID_distance1(nsub) - sub_mean(nsub) + group_mean;
%     Norm_mid2(nsub) = MEAN_MID_distance2(nsub) - sub_mean(nsub) + group_mean;
%     Norm_mid3(nsub) = MEAN_MID_distance3(nsub) - sub_mean(nsub) + group_mean;
%     Norm_mid4(nsub) = MEAN_MID_distance4(nsub) - sub_mean(nsub) + group_mean;
%     Norm_high0(nsub) = MEAN_HIGH_distance0(nsub) - sub_mean(nsub) + group_mean;
%     Norm_high1(nsub) = MEAN_HIGH_distance1(nsub) - sub_mean(nsub) + group_mean;
%     Norm_high2(nsub) = MEAN_HIGH_distance2(nsub) - sub_mean(nsub) + group_mean;
%     Norm_high3(nsub) = MEAN_HIGH_distance3(nsub) - sub_mean(nsub) + group_mean;
%     Norm_high4(nsub) = MEAN_HIGH_distance4(nsub) - sub_mean(nsub) + group_mean;
% end