

lambda_out <- read.table("../sfa_outputs/Deng2014PPS/sqrt_deng_counts_lambda.out")
f_out <- t(read.table("../sfa_outputs/Deng2014PPS/sqrt_deng_counts_F.out"))

indices_mat <- SFA.ExtractTopFeatures(f_out, top_features = 10, options="min")
