# [Evaluating population viability and efficacy of conservation management using integrated population models](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2664.13080)

### Sarah P. Saunders, Francesca J. Cuthbert, and Elise F. Zipkin

### Journal of Applied Ecology

### Please contact the first author for questions about the code or data: Sarah P. Saunders (sarahpsaunders@gmail.com)
__________________________________________________________________________________________________________________________________________
## Abstract
1. Predicting population responses to environmental conditions or management scenarios is a fundamental challenge for conservation. Proper consideration of demographic, environmental, and parameter uncertainties is essential for projecting population trends and optimal conservation strategies. 
2. We developed a coupled integrated population model-Bayesian population viability analysis (IPM-BPVA) to assess the (i) impact of demographic rates (survival, fecundity, immigration) on past population dynamics; (ii) population viability 10 years into the future; and (iii) efficacy of possible management strategies for the federally endangered Great Lakes piping plover (*Charadrius melodus*) population. 
3. Our model synthesizes long-term population survey, nest monitoring, and mark-resight data, while accounting for multiple sources of uncertainty. We incorporated latent abundance of eastern North American merlins (*Falco columbarius*), a primary predator of adult plovers, as a covariate on adult survival via a parallel state-space model, accounting for the influence of an imperfectly observed process (i.e. predation pressure) on population viability. 
4. Mean plover abundance increased from 18 pairs in 1993 to 75 pairs in 2016, but annual population growth (λ_t) was projected to be 0.95 (95% CI: 0.72 – 1.12), suggesting a potential decline to 67 pairs within ten years. Without accounting for an expanding merlin population, we would have concluded that the plover population was projected to increase (λ_t) = 1.02; 95% CI: 0.94 – 1.09) to 91 pairs by 2026. We compared four conservation scenarios: (1) no proposed management; (2) increased control of chick predators (e.g. Corvidae, Laridae, mammals); (3) increased merlin control; and (4) simultaneous chick predator and merlin control. Compared to the null scenario, chick predator control reduced quasi-extinction probability from 11.9% to 8.7%, merlin control more than halved (3.5%) the probability, and simultaneous control reduced quasi-extinction probability to 2.6%. 
5. *Synthesis and applications*. Piping plover recovery actions should consider systematic predator control, rather than current ad hoc protocols, especially given the predicted increase in regional merlin abundance. This approach of combining integrated population models with Bayesian population viability analysis to identify limiting components of the population cycle and evaluate alternative management strategies for conservation decision-making shows great utility for aiding recovery of threatened populations.

## Data
*Please see Dryad for additional data files, http://dx.doi.org/10.5061/dryad.j2906

## Code
1. [PIPL_IPM_Saunders.R](https://github.com/zipkinlab/Saunders_etal_2018_JAE/blob/master/PIPL_IPM_Saunders.R) includes code for integrated population model  
2. [PIPL_IPM_BPVA_Saunders.R](https://github.com/zipkinlab/Saunders_etal_2018_JAE/blob/master/PIPL_IPM_BPVA_Saunders.R) includes code for extending the IPM into a BPVA  
3. [PIPL_IPM_BPVA_Mgmt_Saunders.R](https://github.com/zipkinlab/Saunders_etal_2018_JAE/blob/master/PIPL_IPM_BPVA_Mgmt_Saunders.R) includes code for extending the IPM-BPVA to include evaluation of 4 management scenarios
