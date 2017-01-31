#Strata Module - Provides an object for stratification analysis; accepts array of contingency tables
#Written by Jordan Poles
#Provided under the MIT License (c) 2017 Jordan Poles 

from contingency import *
class Strata:
    def __init__(self, tables):
        self.tables = tables;
    def pool(self):
        a = sum([i.a for i in self.tables])
        b = sum([i.b for i in self.tables])
        c = sum([i.c for i in self.tables])
        d = sum([i.d for i in self.tables])
        return ContTable(a, b, c, d, 0)
    def woolf_correction(self, CI=.95, verbose=1):
        weights = [1/((1/(i.a+.5)) + (1/(i.b+.5)) + (1/(i.c+.5)) + (1/(i.d+.5)))  for i in self.tables]
        adj_logOR = sum([w*i.ss_odds_ratio(verbose=0)["logOR"] for i, w in zip(self.tables, weights)])/sum(weights)
        adj_OR = exp(adj_logOR)
        adj_logORVar = 1/sum(weights)
        z = -norm.interval(CI)[0]
        upper = exp(adj_logOR + z * sqrt(adj_logORVar))
        #Lower bound
        lower = exp(adj_logOR - z * sqrt(adj_logORVar))
        if verbose > 1:
            print("Individual ORs:", [i.odds_ratio(verbose=0)["OR"] for i in tables])
        if verbose:
            print("Woolf adj. OR:", adj_OR)
            if CI:
                print("=> {:d}% CI using logOR: {:f}".format(floor(CI*100), adj_logOR)) 
                print("=> => logOR Var:", adj_logORVar)
                print("=> => Upper Bound:", upper)
                print("=> => Lower Bound:", lower)
        return adj_OR, adj_logOR, (lower, upper)     
    def mantel_haenszel_correction(self, CI=.95, verbose=1):
        weights = [i.b*i.c/i.num_tot for i in self.tables]
        adj_OR = sum([i.a*i.d/i.num_tot for i in self.tables])/sum(weights)
        adj_logOR = log(adj_OR)

        first_num = sum([((i.a + i.d)/i.num_tot)*((i.a*i.d)/i.num_tot) for i in self.tables])
        first_denom = 2*sum([((i.a*i.d)/i.num_tot) for i in self.tables])**2
        second_num = sum([((i.a + i.d)/i.num_tot)*((i.b*i.c)/i.num_tot) + ((i.b + i.c)/i.num_tot)*((i.a*i.d)/i.num_tot) for i in self.tables])
        second_denom = 2*sum([((i.a*i.d)/i.num_tot) for i in self.tables])*sum([((i.b*i.c)/i.num_tot) for i in self.tables])
        third_num = sum([((i.b + i.c)/i.num_tot)*((i.b*i.c)/i.num_tot) for i in self.tables])
        third_denom = 2*sum([((i.b*i.c)/i.num_tot) for i in self.tables])**2
        adj_logORVar = (first_num/first_denom) + (second_num/second_denom) + (third_num/third_denom);

        z = -norm.interval(CI)[0]
        upper = exp(adj_logOR + z * sqrt(adj_logORVar))
        #Lower bound
        lower = exp(adj_logOR - z * sqrt(adj_logORVar))
        if verbose > 1:
            print("Individual ORs:", [i.odds_ratio(verbose=0)["OR"] for i in self.tables])
        if verbose:
            print("Mantel-Haenszel adj. OR:", adj_OR)
            if CI:
                print("=> {:d}% CI using logOR: {:f}".format(floor(CI*100), adj_logOR)) 
                print("=> => logOR Var:", adj_logORVar)
                print("=> => Upper Bound:", upper)
                print("=> => Lower Bound:", lower)
        return adj_OR, adj_logOR, (lower, upper)
    def cmh_test(self, alpha=.05, verbose=1):
        Ai = [(i.a + i.b)*(i.a + i.c)/i.num_tot for i in self.tables]
        Vi = [(i.a + i.b)*(i.c + i.d)*(i.a + i.c)*(i.b + i.d)/(i.num_tot**2)/(i.num_tot-1) for i in self.tables]
        test_statistic = (sum([i.a for i in self.tables]) - sum(Ai))**2/sum(Vi)
        pval = 1-chi2.cdf(test_statistic, 1)
        if verbose:
            print("Cochran-Mantel-Haenszel Test:")
            print("=> test statistic =", test_statistic)
            print("=> p-value ≈", pval)
        if verbose > 1:
            if pval < alpha:
                #Reject null hypothesis
                print("=> Cochran-Mantel-Haenszel Test suggests an association between D & E .")
            else:
                #Stick with null hypothesis
                print("=> Cochran-Mantel-Haenszel Test indicates independence between D & E controling for stratified confounder.")
            print("")
        return pval
