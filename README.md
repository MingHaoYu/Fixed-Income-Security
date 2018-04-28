# Fixed-Income-Security
One factor short-term rates: Vasicek model and CIR model

1. Assume the dynamics of the short-term interest rate, under the risk-neutral measure, are given by the following SDE (Vasicek model):

          𝑑𝑟𝑡=𝜅(𝑟̅−𝑟𝑡)𝑑𝑡+𝜎𝑑𝑊𝑡

          with 𝑟0=5%,𝜎=10%,𝜅=0.82,𝑟̅=5%

(a) Use Monte Carlo Simulation (assume each time step is a day) to find the price of a pure discount bond, with Face Value of $1,000, maturing in 𝑇=0.5 years (at time 𝑡=0)

          𝑃(𝑡,𝑇)=𝔼𝑡∗[$1,000∗𝑒𝑥𝑝(−∫𝑟(𝑠)𝑑𝑠𝑇𝑡)]

(b) Use Monte Carlo Simulation to find the price of a coupon paying bond, with Face Value of $1,000, paying semiannual coupons of $30, maturing in 𝑇=4 years where c={30, 30, 30, 30, 30, 30, 30, 1030}, and T={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}

          𝑃(0,𝐶,𝑇)=𝔼0∗[Σ𝐶𝑖8𝑖=1∗𝑒𝑥𝑝(−∫𝑟(𝑠)𝑑𝑠𝑇𝑖0)]

(c) Use Monte Carlo Simulation to find the price of a European Call option on the pure discount bond in part (a). The option matures in 3 months and has a strike price of 𝐾=$980. Use the explicit formula for the underlying bond price (only for the bond price)

(d) Use Monte Carlo Simulation to find the price of a European Call option on the coupon paying bond in part (b). The option matures in 3 months and has a strike price of 𝐾=$980. Use Monte Carlo simulation for pricing the underlying bond.

(e) Find the price of a European Call option of part (d) by using the explicit formula for the underlying bond price, and reconcile the findings with the ones of part (d).

2. Assume the dynamics of the short-term interest rate, under the risk-neutral measure, are given by the following SDE (CIR model):

          𝑑𝑟𝑡=𝜅(𝑟̅−𝑟𝑡)𝑑𝑡+𝜎√𝑟𝑡𝑑𝑊𝑡

          with 𝑟0=5%,𝜎=12%,𝜅=0.92,𝑟̅=5.5%.
          
(a) Use Monte Carlo Simulation to find at time 𝑡=0 the price 𝑐(𝑡,𝑇,𝑆) of a European Call option, with strike price of 𝐾=$980, maturity of 𝑇=0.5 years on a Pure Discount Bond with Face Value of $1,000, that matures in 𝑆=1 year:

          𝑐(𝑡,𝑇,𝑆)=𝔼𝑡∗[𝑒𝑥𝑝(−∫𝑟(𝑢)𝑑𝑢𝑇𝑡)∗max (𝑃(𝑇,𝑆)−𝐾,0)]
          
(b) Use the Implicit Finite-Difference Method to find at time 𝑡=0 the price 𝑐(𝑡,𝑇,𝑆) of a European Call option, with strike price of 𝐾=$980, maturity of 𝑇=0.5 years on a Pure Discount Bond with Face Value of $1,000, that matures in 𝑆=1 year. The PDE is given as

          𝜕𝑐/𝜕𝑡+1/2𝜎^2𝑟𝜕^2𝑐/𝜕𝑟^2+𝜅(𝑟̅−𝑟)𝜕𝑐/𝜕𝑟−𝑟𝑐=0
          
          with 𝑐(𝑇,𝑇,𝑆)=max(𝑃(𝑇,𝑆)−𝐾,0), and 𝑃(𝑇,𝑆) is computed explicitly.
          
(c) Compute the price 𝑐(𝑡,𝑇,𝑆) of the European Call option above using the explicit formula, and compare it to your findings in parts (a) and (b) and comment on your findings.

3. Assume the dynamics of the short-term interest rate, under the risk-neutral measure, are given by the following system of SDE (G2++ model):

          {𝑑𝑥𝑡=−𝑎𝑥𝑡𝑑𝑡+𝜎𝑑𝑊𝑡1
          {𝑑𝑦𝑡=−𝑏𝑦𝑡𝑑𝑡+𝜂𝑑𝑊𝑡2
          {𝑟𝑡=𝑥𝑡+𝑦𝑡+𝜙𝑡
          
          𝑥0=𝑦0=0, 𝜙0=𝑟0=3%, 𝑑𝑊𝑡1𝑑𝑊𝑡2=𝜌𝑑𝑡,𝜌=0.7, 𝑎=0.1,𝑏=0.3,𝜎=3%,𝜂=8%. Assume 𝜙𝑡=𝑐𝑜𝑛𝑠𝑡=3% for any 𝑡≥0.
          
Use Monte Carlo Simulation to find at time 𝑡=0 the price 𝑝(𝑡,𝑇,𝑆) of a European Put option, with strike price of 𝐾=$950, maturity of 𝑇=0.5 years on a Pure Discount Bond with Face value of $1,000, that matures in 𝑆=1 year. Compare it with the price found by the explicit formula and comment on it.
          
