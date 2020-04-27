
## IBS538 Capstone Project

Astrocytes are the most abundant glial population in the human brain and play an active role in neuronal circuit formation and function through a plethora of functions, including synaptogenesis, calcium dynamics, and neurotransmitter modulation. We have shown that human astrocytes undergo a profound transcriptomic transformation from immature proliferative cells to quiescent mature astrocytes and that this includes the up- and down- regulation of genes closely associated with astrocyte calcium signaling <https://www.cell.com/neuron/fulltext/S0896-6273(15)01019-3>. However, it remains unknown if these transcriptomic changes are accompanied by corresponding functional changes.

We are able to test this hypothesis in a well-controlled setup using hiPSC-derived 3D cortical organoids, an in vitro model system, which we have shown captures the transcriptomic maturation changes that we see in vivo <https://www.cell.com/neuron/fulltext/S0896-6273(17)30683-9>. We will sample astrocytes from 3D cortical organoids at four time points (days in culture) that span the window of maturation and then measure the peak Ca amplitude (dF/F) in response to the addition of 100uM of glutamate (a known stimulant of astrocyte Ca signaling). 

Given that genes closely associated with calcium signlaing change as astrocytes mature, I hypothesize that if astrocytes are sampled throughout critical developmental time points, then the magnitude of the Ca spikes should change over time. The independent variable in this experiment is time in units of days in culture. While time is usually a continuous/measured variable, I will treat it as a discrete factorial variable in this analysis. The dependent variable is peak Ca amplitude (dF/F), which is a continuous/measured variable.

Null hypothesis: Ca amplitude variance associated with time (days in culture) is less than or equal to the residual variance.
Alternate hypothesis: Ca amplitude variance associated with time (days in culture) is greater than the residual variance.

To test this, I will run a one-way ANOVA with related measures because this design has only one predictor variable (day in culture) that is imposed at four levels (100, 150, 200, and 250 days in culture) and because the measurements of my dependent variable (peak amplitude) are intrinsically linked for all four levels of time, within each of four independent replicates. 

Our previous research shows that astrocytes start to form in our organoid system around 80 days in culture and will undergo maturation until 250 days in culture. For this reason, I will serially test astrocytes at four time points: 100, 150, 200, and 250 days in culture, such that measurements of Ca signaling for each independent replicate are intrinsically linked. While our organoids generally exhibit uniform development across various hiPSC lines (derived from different healthy human subjects), I will use four different hiPSC lines (four independent replicates), which was determined through a monte carlo power analysis to give us 99% power for a one-way ANOVA with related measures test. The unit of measurement for the Ca peak amplitude (dF/F) is a commonly used measurement for this type of assay, which uses a fluorescent dye that interacts with intracellular Ca. dF/F details the change in fluorescent intensity between peak and baseline fluorescent intensity to account for any baseline variations. Cells will be imaged for five minutes post-stimulation (addition of glutamate), to enusre enough time for peak Ca signal. This is a standard measurement endpoint for this assay that is well-established in the field.  

```{r message=FALSE, warning=FALSE}
library(tidyverse)
library(viridis)
library(ez)
```

### Simulation data

```{r}
n <- 4
sd <- 0.3
means <- c(0.75, 1.5, 1.8, 2.0)
### end of parameter changes

## you shouldn't need to mess with the rest of this script. you're welcome.

datamaker <- function(n, means, sd) {
  
  params <- data.frame(n, means, sd)
  
  data <- as.tibble(apply(params, 1, function (x) rnorm(x[1], x[2], x[3])))
  
  data <- pivot_longer(data, everything(), 
                       names_to =c("cell_line"),
                       values_to="amplitude") %>% 
    mutate(cell_line=as.factor(paste(
      rep(letters[1:4], each = 4))),
      day=as.factor(paste(
      rep(c("100","150","200", "250"), times = 4))))
  
  data
}

data <- datamaker(n, means, sd)

data
```

```{r}

ggplot(data, aes(as.factor(day), amplitude, color=cell_line, group=cell_line))+
  scale_color_viridis(discrete=T)+
  geom_point(size=3)+
  geom_line() +
  xlab("Day in culture")+
  ylab("Ca peak amplitude (dF/F)")+
  ggtitle("Peak Ca Amplitude in Astrocytes Over Time") +
  theme(plot.title = element_text(size=10, hjust = 0.5))
```

### Monte Carlo analysis

```{r}
sims=100

pval <- replicate(
  sims, {
 
    sample.df <- datamaker(n, means, sd)
    
   sim.ezaov <- ezANOVA(data=sample.df, 
                dv=amplitude, 
                wid=cell_line, 
                within = day,
            detailed=F)
  
  pval <- sim.ezaov$ANOVA[1,5]
    
    }
  )

pwr.pct <- sum(pval<0.05)/sims*100
paste(pwr.pct, sep="", "% power. Change 'n' in your initializer for higher or lower power.")

ggplot(data.frame(pval))+
  geom_histogram(aes(pval), color="#d28e00")+
  labs(x="p-value", title="p-value distribution")
```

**A monte carlo power analysis for related measures indicates that the n=4 used in the above experiment will give us 99-100% power. This level of ANOVA power should be enough to detect an effect of time on peak Ca signal in astrocytes.**



