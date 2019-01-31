[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/stineb/rsofun?branch=master&svg=true)](https://ci.appveyor.com/project/stineb/rsofun)
<a href="https://www.buymeacoffee.com/H2wlgqCLO" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" height="21px" ></a>
[![Github All Releases](https://img.shields.io/github/downloads/atom/atom/total.svg)]()


## Installation

### Development release
To install and load the rpmodel package (development release) run the following command in your R terminal: 
```r
if(!require(devtools)){install.packages(devtools)}
devtools::install_github( "stineb/rpmodel" )
library(rpmodel)
```

### Stable release
rpmodel is not yet available on CRAN. We're working on it.

## Theory

The P-model predicts an optimal ratio of $c_i : c_a$, termed as $\chi$, that balances the costs associated with maintaining the transpiration ($E$) stream and the carboxylation capacity $V_{\text{cmax}}$. It can therefore be used to simulate the acclimation of the photosynthetic machinery to its environment - a mechanism that happens at a time scale of several days to months. At its core, it provides a solution for the optimality criterium
$$
a \; \frac{\partial (E/A)}{\partial \chi} = -b \; \frac{\partial (V_{\mathrm{cmax}}/A)}{\partial \chi}  \;\;\;\;\;\;\;\;\;\;\;\;(1)
$$
The optimal $\chi$ solves the above equation and, with $E = 1.6 g_s D$, $A = g_s (1-\chi)$, and using the Rubisco-limited assimilation rate:
$$
A = A_C = V_{\mathrm{cmax}} \; \frac{\chi\;c_a-\Gamma^{\ast}}{\chi\;c_a + K}
$$ 
is given by:
$$
\chi = \frac{\Gamma^{\ast}}{c_a} + \left(1- \frac{\Gamma^{\ast}}{c_a}\right)\;\frac{\xi}{\xi + \sqrt{D}}
$$
with 
$$
\xi = \sqrt{\frac{b(K+\Gamma^{\ast})}{1.6\;a}}
$$
The unit cost ratio $b/a$ is also referred to as $\beta$. 

A more complete description of the model theory is given in Wang et al. (2017)
and in Stocker et al. (2019). The basic idea was presented in Prentice et al.
(2014)

## Example P-model run

So much for the theory. Let's run the P-model, without $J_{\text{max}}$ limitation, for one set of inputs, being temperature, PPFD, VPD, CO$_2$, elevation, and fAPAR.

To do so, run the `rpmodel()` function from the rsofun package:
```{r, message=FALSE, warning=FALSE}
library(rsofun)
library(dplyr)
# modified seq() function to get a logarithmically spaced sequence
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

## Set parameters
beta <- 146          # unit cost ratio a/b
gamma <- 0.105       # unit cost ratio c/b
kphio <- 0.05        # quantum yield efficiency
c_molmass <- 12.0107 # molar mass, g / mol

## Define environmental conditions
tc <- 20             # temperature, deg C
ppfd <- 300          # mol/m2/d
vpd  <- 1000          # Pa
co2  <- 400          # ppm
elv  <- 0            # m.a.s.l.
fapar <- 1           # fraction  

out_analytical <- rsofun::rpmodel( 
  tc             = tc,
  vpd            = vpd,
  co2            = co2,
  elv            = elv,
  kphio          = kphio,
  beta           = beta,
  fapar          = fapar,
  ppfd           = ppfd,
  method_optci   = "prentice14",
  method_jmaxlim = "none",
  do_ftemp_kphio = FALSE 
  )
```

The function returns a list of variables (see also man page by `?rpmodel`), including $V_{\mathrm{cmax}}$, $g_s$, and all the parameters of the photosynthesis model ($K$, $\Gamma^{\ast}$), which are all internally consistent, as can be verified for...
$$
c_i = c_a - A / g_s = \chi c_a
$$

```{r}
print( out_analytical$ci )
print( out_analytical$ca - (out_analytical$gpp / c_molmass) / out_analytical$gs )
print( out_analytical$ca * out_analytical$chi )
```
Yes. 

And for...
$$
A = V_{\text{cmax}} \frac{c_i-\Gamma^{\ast}}{c_i + K} = \phi_0 I_{\text{abs}} \frac{c_i-\Gamma^{\ast}}{c_i + 2 \Gamma^{\ast}} = g_s (c_a - c_i)
$$

```{r}
print( out_analytical$gpp / c_molmass )
print( out_analytical$vcmax * (out_analytical$ci - out_analytical$gammastar) / (out_analytical$ci + out_analytical$kmm ))
print( out_analytical$gs * (out_analytical$ca - out_analytical$ci) )

print( kphio * ppfd * fapar * (out_analytical$ci - out_analytical$gammastar) / (out_analytical$ci + 2 * out_analytical$gammastar ))
```
Yes.

## Author and contact

Benjamin Stocker
benjamin.stocker@gmail.com

## References

Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancingthe costs of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology, Ecology  Letters,  17,  82–91, 10.1111/ele.12211, http://dx.doi.org/10.1111/ele.12211, 2014.

Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J., and Peng, C.:  Towards a universal model for carbon dioxide uptake byplants, Nat Plants, 3, 734–741, 2017.
