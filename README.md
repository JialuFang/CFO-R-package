# CFO-R-package

## CFO_aCFO

### CFO.next

**Description**

Use the function to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.

**Usage** 

CFO.next(phi, cys, cns, alp.prior=phi, bet.prior=1-phi, cover.doses)

**Parameters**

phi: the target DLT rate

cys: the current number of DLTs observed in patients for the left, current, and right dose levels.

cns: the current number of patients for the left, current, and right dose levels.

alp.prior,bet.prior: the parameters of the prior distribution for the true DLT rate at any dose level.
                    This prior distribution is set to Beta( \code{alpha.prior}, \code{beta.prior}).
                    The default value is \code{phi} and \code{1-phi}.
                    
cover.doses: whether the dose level (left, current and right) is over-toxic or not. 
            The value is set as 1 if the dose level is overly toxicity; otherwise, it is set to 0.
            
**Details**

The CFO design determines the dose level for the next cohort by assessing evidence from the current 
dose level and its adjacent levels. This evaluation is based on odds ratios denoted as \eqn{O_k}, where 
k = L, C, R represents left, current, and right dose levels. Additionally, we define \eqn{\overline{O}_k = 1/O_k}. 
The ratio \eqn{O_C / \overline{O}_{L}} indicates the inclination for de-escalation, while \eqn{\overline{O}_C / O_R} 
quantifies the tendency for escalation. Threshold values \eqn{\gamma_L} and \eqn{\gamma_R} are chosen to 
minimize the probability of making incorrect decisions.The decision process is summarized in Table 1
of Jin and Yin (2022).
An overdose control rule is implemented to ensure patient safety. If the data suggest excessive 
toxicity at the current dose level, we exclude that level and those higher levels. Two scenarios 
lead to a decision on one side only: when the current dose is at the boundary (the first or last dose level) 
or when higher dose levels have been eliminated.

**Values**

The \code{CFO.next()} function returns a list object comprising the following elements: the target DLT 
rate ($target), the current counts of DLTs and patients for the left, current, and right dose levels ($cys and $cns), 
the decision for whether to move to the left or right dose level for the next cohort ($decision), and the 
corresponding index ($index). Specifically, 1 corresponds to de-escalation, 2 corresponds to staying at the 
current dose, and 3 corresponds to escalation.

**Author**

Jialu Fang and Wenliang Wang

**References**

Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.

**Examples**

determine the dose level for the next cohort of new patients

cys <- c(0,1,0); cns <- c(3,6,0)

CFO.next(phi=0.2, cys=cys, cns=cns, alp.prior=0.2, bet.prior=0.8, cover.doses=c(0,0,0))

### aCFO.next

**Description**

In aCFO design, use the function to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.

**Usage** 

aCFO.next(phi, tys, tns, alp.prior, bet.prior, tover.doses, cidx)

**Parameters**

phi: the target DLT rate

tys: the current number of DLTs observed in patients for all dose levels.

tns: the current number of patients for all dose levels.

alp.prior,bet.prior: the parameters of the prior distribution for the true DLT rate at any dose level.
                     This prior distribution is set to Beta( \code{alpha.prior}, \code{beta.prior}). 
                     The default value is \code{phi} and \code{1-phi}.
                     
tover.doses: whether the dose level (from the first to last dose level) is over-toxic or not. 
            The value is set as 1 if the dose level is overly toxicity; otherwise, it is set to 0.
            
cidx: dose level for current cohort
            
**Details**

The aCFO design design incorporate the dose information of all positions (from the lowest to the 
highest dose levels) into the trial decision-making. Prior to assigning dose levels for new patient 
cohorts, aCFO compares the evidence from the current dose level with all doses to its left and right. 
This design is rooted in the odds ratio, specifically \eqn{O_C / \overline{O}_{J}} and 
\eqn{\overline{O}_C / O_R} from the CFO design. By aggregating odds ratios from the left and right sides, 
it forms two collective statistics: \eqn{ {\rm OR}_L =\sum_{i=1}^{J} O_C/ \overline{O}_{L_i} } for dose 
de-escalation (movement to the left) and \eqn{ {\rm OR}_R = \sum_{i=1}^{H} \overline{O}_C / O_{R_i} } 
for dose escalation (movement to the right), where J and H represent the counts of doses on the left and 
right sides of the current dose, respectively. For the new statistic \eqn{ {\rm OR}_L } and \eqn{ {\rm OR}_R }, 
their corresponding thresholds are derived by summing up its individual thresholds \eqn{\gamma_{L_i}} and 
\eqn{\gamma_{R_i}}, i.e., \eqn{\sum_{i=1}^{J}\gamma_{L_i}} and \eqn{\sum_{i=1}^{H}\gamma_{R_i}}. 

Besides，The aCFO design retains the same early stopping criteria as the CFO design. Overall，while preserving 
the nature of the CFO design (model-free and calibration-free), the aCFO designs enhance the efﬁciency by 
incorporating more dose information.

**Values**

The \code{aCFO.next()} function returns a list object comprising the following elements: the target DLT 
rate ($target), the current number of DLTs and patients for the left, current, and right dose levels ($cys and $cns), 
the decision in the aCFO design of whether to move to the left or right dose level for the next cohort ($decision), and the 
corresponding index ($index). Specifically, 1 corresponds to de-escalation, 2 corresponds to staying at the 
current dose, and 3 corresponds to escalation.

**Author**

Jialu Fang 

**References**

Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.

**Examples**

determine the dose level for the next cohort of new patients

tys <- c(0,0,1,0,0,0,0); tns <- c(3,3,6,0,0,0,0)

aCFO.next(phi=0.2, tys=tys, tns=tns, alp.prior=0.2, bet.prior=0.8, tover.doses=c(0,0,0,0,0,0,0), cidx=3)

### aCFO.simu

**Description**

Use this function to find the maximum tolerated dose (MTD) for a single Calibration-Free Odds (CFO) or accumulative CFO (aCFO) trial.

**Usage** 

 aCFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3, alp.prior = phi, bet.prior = 1 - phi, seed=100, accumulation = FALSE)

**Parameters**

phi: the target DLT rate

p.true: the true DLT rates under the different dose levels

ncohort: the total number of cohorts

init.level: the dose level assigned to the first cohort. The default value \code{init.level} is 1.

cohortsize: the sample size in each cohort

alp.prior,bet.prior: the parameters of the prior distribution for the true DLT rate at any dose level.
                    This prior distribution is set to Beta(\code{alpha.prior}, \code{beta.prior}). 
                    The default value is \code{phi} and \code{1-phi}.
                    
seed: the random seed for simulation

accumulation set \code{accumulation=FALSE} to conduct the CFO design; set \code{accumulation=TRUE} to conduct the aCFO design.
            
**Details**

The \code{aCFO.simu()} function is designed to determine the Maximum Tolerated Dose (MTD) for a single CFO or aCFO 
trial. If \code{accumulation = FALSE}, this trial corresponds to the CFO design. If \code{accumulation = TRUE}, it 
corresponds to the aCFO design. 
Given the toxicity outcomes from previous cohorts, each cohort is sequentially assigned to the most suitable dose
level based on the decision rule. Early stopping criteria are incorporated into the CFO design to ensure patient 
safety and benefit. If there is substantial evidence indicating that the current dose level exhibits excessive 
toxicity (\eqn{\Pr(p_C > \phi|x_C, m_C \geq 3) > 0.95}), we exclude the current dose level as well as higher dose 
levels from the trial. Upon the predefined maximum sample size is reached or the lowest dose level is overly toxicity, 
the experiment is concluded, and the MTD is determined using isotonic regression.

**Values**

The \code{aCFO.simu()} function returns a list object comprising the following components: the target DLT 
rate ($target), the actual DLT rates under different dose levels ($p.true), the selected MTD ($MTD), 
the total number of DLTs and patients for all dose levels ($DLT.ns and $dose.ns), and the over-toxicity status 
for all dose levels ($over.doses). Specifically, the value of 1 represents over-toxicity at that dose level, 
while the value of 0 indicates safety at that dose level.

**Author**

Jialu Fang

**References**

Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.

**Examples**

phi <- 0.2; ncohort <- 12; cohortsize <- 3

p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)

find the MTD for a single CFO trial

aCFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3, alp.prior = phi, bet.prior = 1 - phi, seed = 100, accumulation = FALSE)

find the MTD for a single aCFO trial
aCFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3, alp.prior = phi, bet.prior = 1 - phi, seed = 100, accumulation = TRUE)

## CFO_lateonset

### lateonset.next

**Description**

Propose the next dose level in the CFO-type and aCFO-type designs with late-onset toxicities, specifically, including 
Time-to-event CFO (TITE-CFO) design, fractional CFO (fCFO), benchmark CFO, TITE-aCFO design, the f-aCFO design 
and benchmark aCFO design.

**Usage** 

lateonset.next(curDose, phi, tau, impute.method, enter.times, dlt.times, current.t, accumulation, doses, tover.doses, simu, add.args=list(alp.prior=phi, bet.prior=1-phi))

**Parameters**

curDose: the current dose level.

phi: the target DLT rate.

tau: maximal assessment window size

impute.method: the imputing method for handling pending DLT data. \code{impute.method = 'frac'} corresponds to 
                   fractional framework. \code{impute.method = 'TITE'} corresponds to time-to-event framework.
                   \code{impute.method = 'No'} implies no use of any imputing method, corresponding to the 
                   benchmark CFO and benchmark aCFO designs.
                   
enter.times: the time at which each subject existing in the trial enters the trial.

dlt.times: the time to DLT for each subject existing in the trial.  If no DLT occurs for a certain subject, 
           \code{dlt.times} is set to 0.
           
current.t: the current time.

accumulation: set \code{accumulation=FALSE} to conduct the CFO-type design; set \code{accumulation=TRUE} to 
              conduct the aCFO-type design.
              
doses: the dose level for each subject existing in the trial.

tover.doses: whether the dose level (from the first to last dose level) is over-toxic or not. 
             The value is set as 1 if the dose level is overly toxicity; otherwise, it is set to 0.
simu: whether simulation or not, if \code{simu=TRUE}, \code{lateonset.next()} also return \code{tover.doses}.
add.args: additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
          and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
          any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).

**Details**

Late-onset outcomes commonly occur in phase I trials involving targeted agents or immunotherapies. As a 
result, the TITE framework and fractional framework serve as two imputation methods to handle pending data 
related to late-onset outcomes. This approach extends the original designs to integrate time information 
for delayed outcomes, leading to the development of TITE-CFO, fCFO, TITE-aCFO, and f-aCFO designs. \cr
In the TITE framework context, an assumption about the distribution of time to DLT must be pre-specified, 
whereas the fractional framework does not require justification for a specific distribution of the time to 
DLT. Consequently, fCFO and f-aCFO adapt to a more diverse range of scenarios.\cr
The function \code{lateonset.next()} also provides the option to set \code{impute.method = "No"} to execute 
the benchmark CFO and aCFO design. These two methods await complete observation of toxicity outcomes for 
the previous cohorts before determining the next dose assignment. This enhances precision but comes at the 
expense of a prolonged trial duration.

**Values**

The \code{lateonset.next()} function returns the recommended dose level for treating the next cohort of 
patients ($dose). if \code{simu=TRUE}, \code{lateonset.next()} also return a vector indicating whether the 
dose level (from the first to last dose level) is over-toxic or not ($tover.doses).

**Author**

Jialu Fang

**References**

Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
\emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.

Jin, H., & Yin, G. (2023). Time‐to‐event calibration‐free odds design: A robust efficient design for 
phase I trials with late‐onset outcomes. \emph{Pharmaceutical Statistics}.

Yin, G., Zheng, S., & Xu, J. (2013). Fractional dose-finding methods with late-onset toxicity in 
phase I clinical trials. \emph{Journal of Biopharmaceutical Statistics}, 23(4), 856-870.

**Examples**



## 2dCFO

### 2dCFO_next

**Description**

Determine the dose combination for the next cohort based on the toxicity outcomes of the enrolled cohorts.

**Parameters**

target: The target DLT rate

npts: A J*K matrix (J<=K) containing the number of patients treated at each dose combination

ntox: A J*K matrix (J<=K) containing the number of patients experienced dose-limiting toxicity at each dose combination

...

**Values**

(j,k): Recommended dose combination for the next cohort

...


### 2dCFO_sim

**Description**

Conduct simulations and get operating characteristics for the 2dCFO design with customized number of iterations.

**Parameters**

target: The target DLT rate

...






