function [ pvalue ] = get_permutation_pval( real_stat, random_stat )
%real stat = observed data result (i.e. accuracy, Fscore, AUC etc)
%random_stat = vector of randomized accuracies

pvalue = (length(find(random_stat>=real_stat))+1) / (length(random_stat)+1);


end

