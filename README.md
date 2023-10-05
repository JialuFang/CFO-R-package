# CFO-R-package

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






