
BiocManager::install('survcomp')
library(survcomp)

idx <- which(!is.na(daysToDeath))
idx

c <- concordance.index(x=riskScore2[idx], 
                       surv.time=daysToDeath[idx], 
                       surv.event=vitalStatus[idx], 
                       #cl=riskGroup[idx],
                       method="noether")

c$c.index
