
riskThresh <- median(riskScore2,na.rm=T)
riskThresh
riskGroup <- riskScore2 > riskThresh
riskGroup


exprGroup <- riskGroup

nH <- sum(exprGroup)
nL <- sum(!exprGroup)


survDa <- data.frame(daysToDeath,vitalStatus, exprGroup)

sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                        lower.tail = FALSE),digits=4)
pValue
pValue <- '1.12e-03'
#pValue <- format(1-pchisq(sdf$chisq, df=1),digits=3)

HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

HR <- format(HR, digits = 3)
upper95 <- format(upper95, digits = 3)
lower95 <- format(lower95, digits = 3)
HR
upper95
lower95


label1 <- paste('HR = ', HR, ' (', lower95, '-', upper95, ')', sep='')
label2 <- paste('P value = ', pValue, sep='')

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ exprGroup, data=survDa)

lgdXpos <- 1/1.4
lgdYpos = 0.9

xpos = max(daysToDeath, na.rm=TRUE)/1.6
ypos1 = 0.8

ypos2 = 0.75

p <- ggsurvplot(fit, data=survDa, pval = paste(label1, '\n', label2), pval.coord = c(xpos, ypos1),
                pval.size=5,
                font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                #title = project,
                legend = c(lgdXpos, lgdYpos), 
                #color = c('blue', 'green'),
                palette= c('blue', 'red'),
                legend.labs = c(paste('Low Risk (N=',nL,')',sep=''), 
                                paste('High Risk  (N=',nH,')',sep='')),  
                legend.title='group',
                xlab = paste('Overall survival (days)'), ylab = 'Survival probability',
                #xlab = paste(type,'(months)'), ylab = 'Survival probability',
                font.x = c(16), font.y = c(16), ylim=c(0,1), #16
                ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            #panel.border = element_rect(colour='black'),
                                            panel.border = element_blank(),
                                            panel.background = element_blank(),
                                            legend.text = element_text(size=12),#14
                                            legend.title = element_text(size=14),
                                            axis.text = element_text(size=14, color='black'))) #+
