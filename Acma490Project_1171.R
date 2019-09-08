## Acma 490 Project
## Nathan Esau
## 
## Produce the tables and figures used in the Report.

# Load the functions used in the report

## Method #1: Use stocins package

# library(devtools)
# install_github("nathanesau/stocins")
# library(stocins)

## Method #2: Use local scripts

source("insuranceModels.R")
source("interestRateModels.R")

# for creating Age Distribution plot
# see https://raw.githubusercontent.com/mrxiaohe/R_Functions/master/functions/bar

source("bar.R") 

# The .csv files have been converted to R data frame files

load("CdnPop2016.rda")
load("MaleMort91.rda")
load("FemaleMort91.rda")
load("Returns91.rda")

# Load other packages

library(lattice)
library(xtable)
library(forecast)

## Fit OU Model

returnsLast5 <- returns91[(length(returns91[,1])-59):length(returns91[,1]),]
row.names(returnsLast5) <- 1:length(returnsLast5[,1])

plot(returnsLast5$Return - 0.05, type = 'l', ylab = 'Return - 0.05',
     xlab = 'Month')

model = arima(returnsLast5$Return - 0.05, order = c(1, 0, 0),
                method = 'CSS', include.mean = FALSE)

armodel = iratemodel(params = list(phi1 = model$coef, sigma = sqrt(model$sigma2),
                                   delta = 0.05, delta0 = 0.0767),
                     "ar1")

# OU Process

irm = iratemodel.convert("ar1", "ou", armodel, 1/12)
irm$delta = 0.05
irm$delta0 = tail(returnsLast5[,2], 1)

arforecast = forecast(model, h = 120)
arforecast$mean = arforecast$mean + 0.05
arforecast$lower = arforecast$lower + 0.05
arforecast$upper = arforecast$upper + 0.05
arforecast$x = arforecast$x + 0.05

pdf("Report/images/arforecast.pdf")
plot(arforecast, axes = FALSE, xlab = "Month", ylab = "Interest Rate",
     ylim = c(0.015, 0.13), main = "")
box()
axis(1, at = seq(10, 160, 50), labels = c(-50, 0, 50, 100))
axis(2)
abline(h = 0.05, lty = 3)
legend('topright', leg = expression(paste(delta, " = 0.05  ")), lty = 3)
dev.off()

## Mortality

mortmale = mortassumptions(params = list(x = 0, table = "MaleMort91"))
mortfemale = mortassumptions(params = list(x = 0, table = "FemaleMort91"))

maleProb = numeric(101)
femaleProb = numeric(101)

for(k in seq(0, 100, 1))
{
  maleProb[k+1] = kdeferredqx(k, mortmale)
  femaleProb[k+1] = kdeferredqx(k, mortfemale)
}

pdf("Report/images/curvedeaths.pdf")

plot(x = seq(0, 100, 1), y = maleProb * 1e+05, type = 'l', xlab = "Age",
     ylab = "Number of Deaths", ylim = c(0, 4000))

lines(x = seq(0, 100, 1), y = femaleProb * 1e+05, type = 'l', lty = 2)

legend('topleft', leg = c("Male", "Female"),
       lty = c(1, 2))

dev.off()

## Analyze Product for Single Life

create3Dplot <- function(table, type)
{
  g <- expand.grid(y = seq(1,10,1), x = seq(25,75,5))
  g$z <- numeric(nrow(g))

  for(i in 1:nrow(g))
  {
    mort = mortassumptions(list(x = g$x[i], table = table))
    term = insurance(list(n = g$y[i], d = 1, e = 0), "isingle", "endow")
    pure = insurance(list(n = g$y[i], d = 0, e = 1), "isingle", "endow")
    
    ins = term
    
    if(type == "endow") {
      epvterm = z.moment(1, term, mort, irm)
      epvpure = z.moment(1, pure, mort, irm)
      ins = insurance(list(n = g$y[i], d = 1, e = pmin(3*epvterm/epvpure,2)), 
       "isingle", "endow")
    }
    else if(type == "pure")
    {
      ins = pure
    }
                 
    g$z[i] = z.sd(ins, mort, irm) / z.moment(1, ins, mort, irm)
  }
  
  wireframe(z ~ y * x, data = g, drape = TRUE, col = 'black', col.regions = 'white',
            aspect = c(1.0, 0.8), colorkey = FALSE, xlab = "n", ylab = "issue age",
            zlab = "", screen = list(z = 340, x = -70), zoom = 0.925,
            scales = list(arrows = FALSE, col = "black", font = 10, cex = 1.0),
            par.settings = list(regions=list(alpha = 0.3), 
                                axis.line = list(col = "transparent")))
}

pdf("Report/images/termSDPlotMale.pdf")
create3Dplot("MaleMort91", "term")
dev.off()

pdf("Report/images/termSDPlotFemale.pdf")
create3Dplot("FemaleMort91", "term")
dev.off()

pdf("Report/images/pureSDPlotMale.pdf")
create3Dplot("MaleMort91", "pure")
dev.off()

pdf("Report/images/pureSDPlotFemale.pdf")
create3Dplot("FemaleMort91", "pure")
dev.off()

n = seq(1,20,1)
pvvar = numeric(length(n))

for(i in 1:length(n))
{
  # add 1e-14 to avoid negative variance due to rounding
  pvvar[i] = pv.var(n[i], irm) + 1e-14
}

pdf("Report/images/varpv.pdf")
plot(x = n, y = sqrt(pvvar), type = 'l', ylab = "sd", xlab = "n")
dev.off()

pdf("Report/images/endowSDPlotMale.pdf")
create3Dplot("MaleMort91", "endow")
dev.off()

pdf("Report/images/endowSDPlotFemale.pdf")
create3Dplot("FemaleMort91", "endow")
dev.off()

## Create an EPV table for different age

epvTable = data.frame(x = seq(25,75,5), EMaleTerm5 = 0,
                      EMalePure5 = 0, EFemaleTerm5 = 0,
                      EFemalePure5 = 0, EMaleTerm10 = 0,
                      EMalePure10 = 0, EFemaleTerm10 = 0,
                      EFemalePure10 = 0)

for(i in 1:length(epvTable$x))
{
  malemort = mortassumptions(list(x = epvTable$x[i], table = "MaleMort91"))
  femalemort = mortassumptions(list(x = epvTable$x[i], table = "FemaleMort91"))
  
  term = insurance(list(n = 5, d = 1, e = 0), "isingle", "endow")
  pure = insurance(list(n = 5, d = 0, e = 1), "isingle", "endow")
  
  epvTable$EMaleTerm5[i] = z.moment(1, term, malemort, irm)
  epvTable$EFemaleTerm5[i] = z.moment(1, term, femalemort, irm)
  epvTable$EMalePure5[i] = z.moment(1, pure, malemort, irm)
  epvTable$EFemalePure5[i] = z.moment(1, pure, femalemort, irm)
  
  term = insurance(list(n = 10, d = 1, e = 0), "isingle", "endow")
  pure = insurance(list(n = 10, d = 0, e = 1), "isingle", "endow")
  
  epvTable$EMaleTerm10[i] = z.moment(1, term, malemort, irm)
  epvTable$EFemaleTerm10[i] = z.moment(1, term, femalemort, irm)
  epvTable$EMalePure10[i] = z.moment(1, pure, malemort, irm)
  epvTable$EFemalePure10[i] = z.moment(1, pure, femalemort, irm)
}

epvTableLatex = xtable(epvTable, digits = c(0,0,4,4,4,4,4,4,4,4),
                       caption = "$EPV$ for term and pure endowment contracts used in pricing",
                       align = "l|l|rrrr|rrrr|")

names(epvTableLatex) = c("Age", "$E$[Term]", "$E$[Pure]", "$E$[Term]", "$E$[Pure]",
                         "$E$[Term]", "$E$[Pure]", "$E$[Term]", "$E$[Pure]")

addtorow <- list(pos = list(-1, -1))
addtorow$command <- c("\\hline\n& \\multicolumn{4}{c|}{5-year} & \\multicolumn{4}{c|}{10-year} \\\\\n",
                      "& \\multicolumn{2}{c}{Male} & \\multicolumn{2}{c|}{Female} & \\multicolumn{2}{c}{Male} & \\multicolumn{2}{c|}{Female} \\\\\n")

fileConn = file("Report/tables/epvTable.tex")
writeLines(print(epvTableLatex,
                 table.placement = "!htpb",
                 include.rownames = FALSE,
                 sanitize.text = function(x){x}, size = "small",
                 add.to.row = addtorow), fileConn)
close(fileConn)

## Create a face amount table for different ages

faceAmountTable = data.frame(x = seq(25,75,5), bmale5 = 1, emale5 = 0, NSPmale5 = 0, bfemale5 = 1,
                             efemale5 = 0, NSPfemale5 = 0, bmale10 = 1, emale10 = 0, NSPmale10 = 0,
                             bfemale10 = 1, efemale10 = 0, NSPfemale10 = 0)

for(i in 1:length(epvTable$x))
{
  faceAmountTable$emale5[i] = min(3 * epvTable$EMaleTerm5[i] / epvTable$EMalePure5[i], 2)
  faceAmountTable$efemale5[i] = min(3 * epvTable$EFemaleTerm5[i] / epvTable$EFemalePure5[i], 2)
  faceAmountTable$emale10[i] = min(3 * epvTable$EMaleTerm10[i] / epvTable$EMalePure10[i], 2)
  faceAmountTable$efemale10[i] = min(3 * epvTable$EFemaleTerm10[i] / epvTable$EFemalePure10[i], 2)
  
  faceAmountTable$NSPmale5[i] = 1 * epvTable$EMaleTerm5[i] + faceAmountTable$emale5[i] * epvTable$EMalePure5[i]
  faceAmountTable$NSPfemale5[i] = 1 * epvTable$EFemaleTerm5[i] + faceAmountTable$efemale5[i] * epvTable$EFemalePure5[i]
  faceAmountTable$NSPmale10[i] = 1 * epvTable$EMaleTerm10[i] + faceAmountTable$emale10[i] * epvTable$EMalePure10[i]
  faceAmountTable$NSPfemale10[i] = 1 * epvTable$EFemaleTerm10[i] + faceAmountTable$efemale10[i] * epvTable$EFemalePure10[i]
  
  # profit and expenses loading
  faceAmountTable$NSPmale5[i] = 1.05 * faceAmountTable$NSPmale5[i]
  faceAmountTable$NSPfemale5[i] = 1.05 * faceAmountTable$NSPfemale5[i]
  faceAmountTable$NSPmale10[i] = 1.05 * faceAmountTable$NSPmale10[i]
  faceAmountTable$NSPfemale10[i] = 1.05 * faceAmountTable$NSPfemale10[i]
}

faceAmountTable[,2:13] = faceAmountTable[,2:13] * 1000

faceAmountTableLatex = xtable(faceAmountTable, digits = c(0,0,0,1,1,0,1,1,0,1,1,0,1,1),
                              caption = paste("Net Single Premium in Thousands of Dollars for Death Benefit of",
                                              "\\$1,000,000 and Survival Benefit satisfying the regulation.",
                                              "A 5\\% profit and expense loading is included in the premium."), 
                              align = "l|l|rrrrrr|rrrrrr|",
                              label = "tab:faceAmountTable")

addtorow <- list(pos = list(-1, -1))
addtorow$command <- c("\\hline\n& \\multicolumn{6}{c|}{5-year} & \\multicolumn{6}{c|}{10-year} \\\\\n",
                      "& \\multicolumn{3}{c}{Male} & \\multicolumn{3}{c|}{Female} & \\multicolumn{3}{c}{Male} & \\multicolumn{3}{c|}{Female} \\\\\n")

names(faceAmountTableLatex) = c("Age", "$d$", "$e$", "NSP", "$d$",
                                "$e$", "NSP", "$d$", "$e$", "NSP", "$d$",
                                "$e$", "NSP")

fileConn = file("Report/tables/faceAmountTable.tex")
writeLines(print(faceAmountTableLatex,
                 include.rownames = FALSE,
                 sanitize.text = function(x){x}, size = "small",
                 add.to.row = addtorow), fileConn)
close(fileConn)

## Investment Return for survival benefit

returnTable = data.frame(x = seq(25,75,5),
                         return5male = log(faceAmountTable$emale5 / 1.05 / faceAmountTable$NSPmale5)/5,
                         return5female = log(faceAmountTable$efemale5 / 1.05 / faceAmountTable$NSPfemale5)/5,
                         return10male = log(faceAmountTable$emale10 / 1.05 / faceAmountTable$NSPmale10)/10,
                         return10female = log(faceAmountTable$efemale10 / 1.05 / faceAmountTable$NSPfemale10)/10)

returnTableLatex = returnTable
returnTableLatex[,2:5] = returnTableLatex[,2:5] * 100
returnTableLatex = xtable(returnTableLatex, caption = "Investment return on the survival benefit conditional on the policyholder surviving the term of the contract.",
                          digits = c(0,0,4,4,4,4), align = "l|l|rr|rr|",
                          label = "tab:returnTable")

names(returnTableLatex) = c("Age", "Return (\\%)", "Return (\\%)", "Return (\\%)", "Return (\\%)")

addtorow <- list(pos = list(-1, -1))
addtorow$command <- c("\\hline\n& \\multicolumn{2}{c|}{5-year} & \\multicolumn{2}{c|}{10-year} \\\\\n",
                      "& \\multicolumn{1}{c}{Male} & \\multicolumn{1}{c|}{Female} & \\multicolumn{1}{c}{Male} & \\multicolumn{1}{c|}{Female} \\\\\n")

fileConn = file("Report/tables/returnTable.tex")
writeLines(print(returnTableLatex,
                 include.rownames = FALSE,
                 sanitize.text = function(x){x}, size = "small",
                 add.to.row = addtorow), fileConn)
close(fileConn)

## Risk for a portfolio of contracts

createSDplot <- function(n, table)
{
  term = insurance(list(n = n, d = 1, e = 0), "isingle", "endow")
  pure = insurance(list(n = n, d = 0, e = 1), "isingle", "endow")
  
  x = seq(25,75,5)
  
  sd1 = numeric(length(x))
  sd10 = numeric(length(x))
  sd100 = numeric(length(x))
  sd1000 = numeric(length(x))
  sdInf = numeric(length(x))
    
  for(i in 1:length(x))
  {
    mort = mortassumptions(list(x = x[i], table = table))
    
    e = 1000 * min(3 * z.moment(1, term, mort, irm) / z.moment(1, pure, mort, irm) , 2)
    d = 1000
    
    endow = insurance(list(n = n, d = d, e = e), "isingle", "endow")
    port1 = insurance(list(single = endow, c = 1), "iport", "endow")
    port10 = insurance(list(single = endow, c = 10), "iport", "endow")
    port100 = insurance(list(single = endow, c = 100), "iport", "endow")
    port1000 = insurance(list(single = endow, c = 1000), "iport", "endow")
    portInf = insurance(list(single = endow, c = 1e18), "iport", "endow")
    
    # insurance risk
    sd1[i] = z.sd(port1, mort, irm) / port1$c
    sd10[i] = z.sd(port10, mort, irm) / port10$c
    sd100[i] = z.sd(port100, mort, irm) / port100$c
    sd1000[i] = z.sd(port1000, mort, irm) / port1000$c
    sdInf[i] = z.sd(portInf, mort, irm) / portInf$c
    
    sd1[i] = sd1[i] / z.moment(1, endow, mort, irm)
    sd10[i] = sd10[i] / z.moment(1, endow, mort, irm)
    sd100[i] = sd100[i] / z.moment(1, endow, mort, irm)
    sd1000[i] = sd1000[i] / z.moment(1, endow, mort, irm)
    sdInf[i] = sdInf[i] / z.moment(1, endow, mort, irm)
  }
  
  plot(x = x, y = sd1, type = 'l', ylim = c(0, max(sd1) * 1.4),
       ylab = "sd/NSP", xlab = "Age")
  lines(x = x, y = sd10, type = 'l', lty = 2)
  lines(x = x, y = sd100, type = 'l', lty = 3)
  lines(x = x, y = sd1000, type = 'l', lty = 4)
  lines(x = x, y = sdInf, type = 'l', lty = 5)
  
  legend('topleft', leg = c("c = 1", "c = 10", "c = 100", "c = 1000", "c = Inf"),
         lty = c(1,2,3,4,5), cex = 0.9, ncol = 5)
  
  return(sdInf)
}

pdf("Report/images/portfolioRiskMale5.pdf")
invriskMale5 = createSDplot(5, "MaleMort91")
names(invriskMale5) = seq(25,75,5)
dev.off()

pdf("Report/images/portfolioRiskFemale5.pdf")
invriskFemale5 = createSDplot(5, "FemaleMort91")
names(invriskFemale5) = seq(25,75,5)
dev.off()

pdf("Report/images/portfolioRiskMale10.pdf")
invriskMale10 = createSDplot(10, "MaleMort91")
names(invriskMale10) = seq(25,75,5)
dev.off()

pdf("Report/images/portfolioRiskFemale10.pdf")
invriskFemale10 = createSDplot(10, "FemaleMort91")
names(invriskFemale10) = seq(25,75,5)
dev.off()

cat("5 year male contract: sd[Z(c)]/E[Z(c)]\n")
print(invriskMale5)

cat("5 year female contract: sd[Z(c)]/E[Z(c)]\n")
print(invriskFemale5)

cat("10 year female contract: sd[Z(c)]/E[Z(c)]\n")
print(invriskMale10)

cat("10 year female contract: sd[Z(c)]/E[Z(c)]\n")
print(invriskFemale10)

## Age Distribution

PopulationTable <- CdnPop2016[6:15,]
PopulationTable$MaleProp <- PopulationTable$MaleProp / sum(PopulationTable$MaleProp)
PopulationTable$FemaleProp <- PopulationTable$FemaleProp / sum(PopulationTable$FemaleProp)
names(PopulationTable) = c("AgeGroup", "Total", "Total", "Prop", "Prop")

MalePopulationTable = PopulationTable[,c(1,2,4)]
MalePopulationTable$Gender = "Male"
FemalePopulationTable = PopulationTable[,c(1,3,5)]
FemalePopulationTable$Gender = "Female"

PopulationTable = rbind(MalePopulationTable, FemalePopulationTable)
PopulationTable$AgeGroup = as.character(PopulationTable$AgeGroup)
PopulationTable$AgeGroup = c("25-29","29-34","35-39","39-44","45-49","50-54",
                              "55-59","60-64","65-69","70-74")
row.names(PopulationTable) = 1:20

PopulationTableAssumed = PopulationTable
PopulationTableAssumed$Prop = rep(c(0.05, 0.075, 0.100, 0.125, 0.1500, 0.1500, 
                                    0.125, 0.100, 0.075, 0.05), 2)
PopulationTableAssumed$Gender = c(rep("Assumed Female", 10),
                                  rep("Assumed Male", 10))

PopulationTable = rbind(PopulationTable, PopulationTableAssumed)

pdf("Report/images/ageDistribution.pdf")

bar(dv = Prop, factors = c(Gender, AgeGroup), dataframe = PopulationTable, 
    errbar = FALSE, ylim = c(0, 0.2), cex.names = 0.8, cex.axis = 0.8, ylab = "Proportion",
    xlab = "Age Group", density = c(30,30,10,10), col = c('red','blue','red','blue'),
    args.legend = list(ncol = 2))

dev.off()

## Risk for a group of contracts

# if getEach is TRUE then get risks for each group individually
getRisk <- function(numPolicies, getEach = FALSE)
{
  PropMale = 0.5 # proportion of endowment contracts issued to males
  PropFemale = 1 - PropMale # proportion of endowment contracts issued to females
  Prop5Year = 0.5 # proportion of endowment contracts with term 5 years
  Prop10Year = 1 - Prop5Year # proportion of endowment contracts with term 10 years
  Prop500 = 0.8 # proportion of endowment contracts with death benefit 500,000
  Prop5000 = 0.2 # proportion of endowment contracts with death benefit 5,000,000

  GroupTable = data.frame(
    Term = c(5, 5, 5, 5, 10, 10, 10, 10),
    DeathBenefit = c(500, 500, 5000, 5000, 500, 500, 5000, 5000),
    Gender = c("Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female"),
    Ages = c("25 to 75", "25 to 75", "25 to 75", "25 to 75", "25 to 75", "25 to 75", "25 to 75", "25 to 75"),
    NumberPolicies = c(numPolicies * PropMale * Prop5Year * Prop500,
                       numPolicies * PropFemale * Prop5Year * Prop500,
                       numPolicies * PropMale * Prop5Year * Prop5000,
                       numPolicies * PropFemale * Prop5Year * Prop5000,
                       numPolicies * PropMale * Prop10Year * Prop500,
                       numPolicies * PropFemale * Prop10Year * Prop500,
                       numPolicies * PropMale * Prop10Year * Prop5000,
                       numPolicies * PropFemale * Prop10Year * Prop5000)
  )

  counts = c(PopulationTableAssumed$Prop[1:10] * GroupTable$NumberPolicies[1],
             PopulationTableAssumed$Prop[11:20] * GroupTable$NumberPolicies[2],
             PopulationTableAssumed$Prop[1:10] * GroupTable$NumberPolicies[3],
             PopulationTableAssumed$Prop[11:20] * GroupTable$NumberPolicies[4],
             PopulationTableAssumed$Prop[1:10] * GroupTable$NumberPolicies[5],
             PopulationTableAssumed$Prop[11:20] * GroupTable$NumberPolicies[6],
             PopulationTableAssumed$Prop[1:10] * GroupTable$NumberPolicies[7],
             PopulationTableAssumed$Prop[11:20] * GroupTable$NumberPolicies[8])

  groups = vector("list", 80)
  mortgroups = vector("list", 80)
  groupNumber = 0
  
  riskTable = data.frame(
    Group = 1:80,
    NumberPH = counts,
    Age = rep(c(27, 32, 37, 42, 47, 52, 57, 62, 67, 72), 8),
    Gender = rep(c(rep("Male", 10), rep("Female", 10)), 4),
    Term = c(rep(5, 40), rep(10, 40)),
    DeathBenefit = rep(c(rep(500,20), rep(5000,20)), 2),
    NPV = numeric(80),
    InsRisk = numeric(80),
    InvRisk = numeric(80),
    TotalSD = numeric(80),
    Ratio = numeric(80)
  ) # only used if getEach is TRUE

  for(term in c(5,10))
  {
    for(deathBenefit in c(500, 5000))
    {
      for(table in c("MaleMort91", "FemaleMort91"))
      {
        for(age in c(27, 32, 37, 42, 47, 52, 57, 62, 67, 72))
        {
          groupNumber = groupNumber + 1
        
          mort = mortassumptions(list(x = age, table = table))
          termins = insurance(list(n = term, d = 1, e = 0), "isingle", "endow")
          pureins = insurance(list(n = term, d = 0, e = 1), "isingle", "endow")
          e = 3 * z.moment(1, termins, mort, irm) / z.moment(1, pureins, mort, irm)
        
          endow = insurance(list(n = term, d = deathBenefit, e = e * deathBenefit), "isingle", "endow")
          endowport = insurance(list(single = endow, c = counts[groupNumber]), "iport", "endow")
        
          groups[[groupNumber]] = endowport
          mortgroups[[groupNumber]] = mort
          
          if(getEach == TRUE)
          {
            currentgroup = insurance(list(endowport), "igroup")
            
            npv = 1.05 * z.moment(1, currentgroup, list(mort), irm) / currentgroup$c
            insrisk = z.insrisk(currentgroup, list(mort), irm) / currentgroup$c^2
            invrisk = z.invrisk(currentgroup, list(mort), irm) / currentgroup$c^2
            totalsd = (insrisk + invrisk)^0.5
            
            riskTable$NPV[groupNumber] = npv
            riskTable$InsRisk[groupNumber] = insrisk
            riskTable$InvRisk[groupNumber] = invrisk
            riskTable$TotalSD[groupNumber] = totalsd
            riskTable$Ratio[groupNumber] = totalsd / npv
          }
        }
      }
    }
  }
  
  if(getEach == FALSE)
  {
    groupins = insurance(groups, "igroup")
    npv = 1.05 * z.moment(1, groupins, mortgroups, irm) / groupins$c
    insrisk = z.insrisk(groupins, mortgroups, irm) / groupins$c^2
    invrisk = z.invrisk(groupins, mortgroups, irm) / groupins$c^2
    totalsd = (insrisk + invrisk)^0.5
  
    return(list(numPolicies = numPolicies, npv = npv, insrisk = insrisk, invrisk = invrisk,
       totalsd = totalsd, ratio = totalsd/npv))
  }
  else
  {
    return(riskTable)
  }
}

groupsTable = getRisk(2000, TRUE)

groupsTableLatex = cbind(groupsTable[1:40,1:6], groupsTable[41:80,1:6])
groupsTableLatex = xtable(groupsTableLatex,
  digits = c(0,0,0,0,0,0,0,0,0,0,0,0,0),
  align = "llrrrrr|lrrrrr",
  caption = "The illustrative portfolio used. The endowment benefit for each group is calculated using Equation 1.",
  label = "tab:groupstable")
names(groupsTableLatex) = rep(c("Group", "Size", "Age", "Gender", "$n$", "$d$"), 2)

fileConn = file("Report/tables/groupsTable.tex")
writeLines(print(groupsTableLatex,
                 include.rownames = FALSE,
                 table.placement = "!htpb",
                 size = "small",
                 sanitize.text = function(x){x}), fileConn)
close(fileConn)

riskTable = data.frame(
  NumPolicies = c(200,1000,2000,4000,20000,1e+18),
  NPV = numeric(6),
  InsRisk = numeric(6),
  InvRisk = numeric(6),
  TotalSD = numeric(6),
  Ratio = numeric(6)
)

for(i in 1:length(riskTable$NumPolicies))
{
  risk = getRisk(riskTable$NumPolicies[i], FALSE)
  riskTable$NumPolicies[i] = risk$numPolicies
  riskTable$NPV[i] = risk$npv
  riskTable$InsRisk[i] = risk$insrisk
  riskTable$InvRisk[i] = risk$invrisk
  riskTable$TotalSD[i] = risk$totalsd
  riskTable$Ratio[i] = risk$ratio
}

riskTable$NumPolicies[6] = "$\\infty$"

riskTableLatex = xtable(riskTable, digits = c(0,0,2,2,2,2,4),
                        caption = "Investment and Insurance Risk for the Illustrative Portfolio",
                        label = "tab:risktable")
names(riskTableLatex) = c("Num Policies", "NSP", "Ins Risk", "Inv Risk", "SD", "SD/NSP")

fileConn = file("Report/tables/riskTable.tex")
writeLines(print(riskTableLatex,
                 include.rownames = FALSE,
                 sanitize.text = function(x){x}, size = "small"),
                 fileConn)
close(fileConn)

## Recalculate risk table using covariance equivalent AR(1) process

model = arima(returnsLast5$Return - 0.05, order = c(1, 0, 0),
              method = 'CSS', include.mean = FALSE)

armodel = iratemodel(params = list(phi1 = model$coef, sigma = sqrt(model$sigma2),
                                   delta = 0.05, delta0 = 0.0767),
                     "ar1")

armodel$phi1 = armodel$phi1^12
armodel$sigma = armodel$sigma * sqrt(12)

irm = armodel

riskTable = data.frame(
  NumPolicies = c(200,1000,2000,4000,20000,1e+18),
  NPV = numeric(6),
  InsRisk = numeric(6),
  InvRisk = numeric(6),
  TotalSD = numeric(6),
  Ratio = numeric(6)
)

for(i in 1:length(riskTable$NumPolicies))
{
  risk = getRisk(riskTable$NumPolicies[i], FALSE)
  riskTable$NumPolicies[i] = risk$numPolicies
  riskTable$NPV[i] = risk$npv
  riskTable$InsRisk[i] = risk$insrisk
  riskTable$InvRisk[i] = risk$invrisk
  riskTable$TotalSD[i] = risk$totalsd
  riskTable$Ratio[i] = risk$ratio
}

riskTable$NumPolicies[6] = "$\\infty$"

riskTableLatex = xtable(riskTable, digits = c(0,0,2,2,2,2,4),
                        caption = "Investment and Insurance Risk for the Illustrative Portfolio",
                        label = "tab:risktableAR1")
names(riskTableLatex) = c("Num Policies", "NSP", "Ins Risk", "Inv Risk", "SD", "SD/NSP")

fileConn = file("Report/tables/riskTableAR1.tex")
writeLines(print(riskTableLatex,
                 include.rownames = FALSE,
                 table.placement = "!htpb",
                 sanitize.text = function(x){x}, size = "small"),
           fileConn)
close(fileConn)

## Recalculate risk table using deterministic interest rate

irm = iratemodel(list(delta = 0.05), "determ")

riskTable = data.frame(
  NumPolicies = c(200,1000,2000,4000,20000,1e+18),
  NPV = numeric(6),
  InsRisk = numeric(6),
  InvRisk = numeric(6),
  TotalSD = numeric(6),
  Ratio = numeric(6)
)

for(i in 1:length(riskTable$NumPolicies))
{
  risk = getRisk(riskTable$NumPolicies[i], FALSE)
  riskTable$NumPolicies[i] = risk$numPolicies
  riskTable$NPV[i] = risk$npv
  riskTable$InsRisk[i] = risk$insrisk
  riskTable$InvRisk[i] = risk$invrisk
  riskTable$TotalSD[i] = risk$totalsd
  riskTable$Ratio[i] = risk$ratio
}

riskTable$NumPolicies[6] = "$\\infty$"

riskTableLatex = xtable(riskTable, digits = c(0,0,2,2,2,2,4),
                        caption = "Investment and Insurance Risk for the Illustrative Portfolio",
                        label = "tab:risktableDeterm")
names(riskTableLatex) = c("Num Policies", "NSP", "Ins Risk", "Inv Risk", "SD", "SD/NSP")

fileConn = file("Report/tables/riskTableDeterm.tex")
writeLines(print(riskTableLatex,
                 include.rownames = FALSE,
                 table.placement = "!htpb",
                 sanitize.text = function(x){x}, size = "small"),
           fileConn)
close(fileConn)

## Sensitivity Test for OU Process

model = arima(returns91$Return - mean(returns91$Return), order = c(1, 0, 0),
              method = 'CSS', include.mean = FALSE)

armodel = iratemodel(params = list(phi1 = model$coef, sigma = sqrt(model$sigma2)),
                     "ar1")

# OU process for 1991 to 2017

irm = iratemodel.convert("ar1", "ou", armodel, 1/12)
irm$delta = mean(returns91$Return)
irm$delta0 = tail(returnsLast5[,2], 1)

arforecast = forecast(model, h = 120)
arforecast$mean = arforecast$mean + mean(returns91$Return)
arforecast$lower = arforecast$lower + mean(returns91$Return)
arforecast$upper = arforecast$upper + mean(returns91$Return)
arforecast$x = arforecast$x + mean(returns91$Return)

pdf("Report/images/arforecast1991.pdf")
plot(arforecast, axes = FALSE, xlab = "Month", ylab = "Interest Rate",
     ylim = c(0.0, 0.25), main = "")
box()
axis(1, at = seq(15, 415, 100), labels = c(-300, -200, -100, 0, 100))
axis(2)
abline(h = mean(returns91$Return), lty = 3)
legend('topright', leg = expression(paste(delta, " = 0.1159 ")), lty = 3)
dev.off()

riskTable = data.frame(
  NumPolicies = c(200,1000,2000,4000,20000,1e+18),
  NPV = numeric(6),
  InsRisk = numeric(6),
  InvRisk = numeric(6),
  TotalSD = numeric(6),
  Ratio = numeric(6)
)

for(i in 1:length(riskTable$NumPolicies))
{
  risk = getRisk(riskTable$NumPolicies[i], FALSE)
  riskTable$NumPolicies[i] = risk$numPolicies
  riskTable$NPV[i] = risk$npv
  riskTable$InsRisk[i] = risk$insrisk
  riskTable$InvRisk[i] = risk$invrisk
  riskTable$TotalSD[i] = risk$totalsd
  riskTable$Ratio[i] = risk$ratio
}

riskTable$NumPolicies[6] = "$\\infty$"

riskTableLatex = xtable(riskTable, digits = c(0,0,2,2,2,2,4),
                        caption = "Investment and Insurance Risk for the Illustrative Portfolio",
                        label = "tab:risktable1991")
names(riskTableLatex) = c("Num Policies", "NSP", "Ins Risk", "Inv Risk", "SD", "SD/NSP")

fileConn = file("Report/tables/riskTable1991.tex")
writeLines(print(riskTableLatex,
                 include.rownames = FALSE,
                 table.placement = "!htpb",
                 sanitize.text = function(x){x}, size = "small"),
           fileConn)
close(fileConn)
