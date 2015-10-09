This folder contains all the existing sit(social influence task)-related models. Below is the model index, in a chronological order based on development history.

RevLearn_

 - RL: simple RL model with Rescorla Wagner update rule.
 
 - RLnc: ficticious RL model, also model the non-chosen choice (-reward - value[non-chosen]).
 - RLnc_2lr: add an additional learning rate to the non-chosen update
 - RLnc_cfa: add an additional conterfactual attention parameter, (-cfa*reward)
 - RLnc_2lr_cfa: putting 2lr and cfa together
 
 - RLcoh: reweight the values according to decision type (with vs. against)
 - RLcoh_2lr
 - RLcoh_cfa
 - RLcoh_2lr_cfa
 
 - RLcumrew: reweight the values according to coplayers' reward history (cumulative reward)
 - RLcumrew_2lr
 - RLcumrew_cfa
 - RLcumrew_2lr_cfa
 
 - RLcoh_modvalue: the coherence reweight goes into the value directly, the same as RLcoh
 - RLcoh_modprob_tempin: the coherence reweight only goes into the softmax transformation, temperature*value + reweight
 - RLcoh_modprob_tempout: the coherence reweight only goes into the softmax transformation, temperature*(value + reweight)
 
 - RLbeta_alt1: weighted sum of otherRewards according to subjects' preference towards co-players
 
 - RLbeta_alt2: use beta_cdf to quantify others' choice tendency, and then used to model otherValue
 - RLbeta_alt2_c_v1: with evidence weight
 - RLbeta_alt2_c_v2: without evidence weight
 
 - RLbeta_alt3: treat each co-player as a separate reinforcer, use a simple RL update to quantify otherValue
 - RLbeta_alt3_p1_v1: step1, get lr and tau for the others, single lr and tau
 - RLbeta_alt3_p2_v1: 1 lr for chosen and non-chosen choices, with c_rep 
 - RLbeta_alt3_p1_v2: step1, get lr and tau for the others, 4 lr and tau
 - RLbeta_alt3_p2_v2: 1 lr for chosen and non-chosen choices, without c_rep; essentially, this should be the same as _p2_v1
 
 - RLbeta_alt4: (preference) weighted sum of cumulative otherRewards to represent otherValue
 - RevLearn_RLbeta_alt4_c_w_v1_1lr: only uses current preference weights and current with/against
 - RevLearn_RLbeta_alt4_c_w_v2_1lr: preference weight sum up to one --> [3 2 1 1]/7
 - RevLearn_RLbeta_alt4_c_w_v3_1lr: v3_without weight on otherValue when update
 - RevLearn_RLbeta_alt4_c_w_v4_1lr: use actual weight and with/against
 - RevLearn_RLbeta_alt4_c_w_v5_1lr: beta1-4 separately model myValue and otherValue (chosen vs. non-chosen)
 - RevLearn_RLbeta_alt4_c_w_v6_1lr: actual weight but only current with/against
 - RevLearn_RLbeta_alt4_c_w_v7_1lr: based on v1, windows size = 4
 - RevLearn_RLbeta_alt4_c_w_v8_1lr: same beta for myV, separate beta for otherV
 - RevLearn_RLbeta_alt4_c_w_v9_1lr: based on v1, windows size = 5
 - RevLearn_RLbeta_alt4_c_w_v10_1lr: based on v1, adding a cfa parameter
 
 
 - _w: weighted coherence, i.e. weight .* with, weight .* against
 - _n: normalised other Value, i.e. othV(c2) / (othV(c2), othV(~c2))
 
 
 