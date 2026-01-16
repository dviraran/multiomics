#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

//' Exact mHG p-value calculation using Dynamic Programming
//'
//' Calculates the exact p-value for the Minimum HyperGeometric statistic
//' using the O(NB) algorithm described in Eden et al. (2007).
//'
//' @param limits IntegerVector of length N, where limits[i] is the maximum
//' allowed number of hits (b) in the top i+1 ranks to remain non-significant.
//' @param N Total number of genes.
//' @param B Number of genes in the gene set.
//' @return The exact p-value.
//' @export
// [[Rcpp::export]]
double mhg_exact_pvalue(IntegerVector limits, int N, int B) {
    if (B == 0) return 1.0;
    if (B > N) stop("B cannot be larger than N");
    if (limits.size() != N) stop("Limits vector size must match N");

    // dp[k] stores the probability of having exactly k hits (1s) at the current step
    // in a random permutation.
    // We use long double to minimize underflow issues during intermediate steps,
    // though for very small p-values the result might still underflow double.
    std::vector<long double> dp(B + 1, 0.0);
    dp[0] = 1.0;

    for (int i = 1; i <= N; ++i) {
        std::vector<long double> next_dp(B + 1, 0.0);
        int max_k = std::min(i, B); 
        int prev_steps = i - 1;
        int remaining_total = N - prev_steps;
        
        // This should not happen if loop logic is correct
        if (remaining_total <= 0 && i < N) break;

        for (int k = 0; k <= max_k; ++k) {
             // Calculate contribution from state (i-1, k) -> choosing 0
             // Transition: next is 0.
             // Prob = (Zeros Left) / (Total Left)
             // Zeros Left = (N - B) - (prev_steps - k)
             int zeros_left = (N - B) - (prev_steps - k);
             
             if (zeros_left > 0) {
                 long double p_zero = (long double)zeros_left / remaining_total;
                 next_dp[k] += dp[k] * p_zero;
             }

             // Calculate contribution from state (i-1, k-1) -> choosing 1
             // Transition: next is 1.
             // Prob = (Ones Left) / (Total Left)
             if (k > 0) {
                 int ones_left = B - (k - 1);
                 if (ones_left > 0) {
                     long double p_one = (long double)ones_left / remaining_total;
                     next_dp[k] += dp[k-1] * p_one;
                 }
             }
        }

        // Apply boundary constraint
        // limits[i-1] is the max allowed hits at step i to NOT have a significant p-value so far.
        // We prune paths that exceed this limit.
        int limit = limits[i-1];
        for (int k = 0; k <= max_k; ++k) {
            if (k > limit) {
                next_dp[k] = 0.0;
            }
        }
        
        dp = next_dp;
    }

    // dp[B] contains the probability that a random permutation is "valid"
    // (i.e., never crosses the significance threshold).
    // The p-value is the probability of the complement (crossing the threshold).
    double p_val = 1.0 - (double)dp[B];
    
    // Ensure numerical range
    if (p_val < 0) p_val = 0.0;
    if (p_val > 1) p_val = 1.0;

    return p_val;
}
