# Schisto IBM   

## Matthew Malishev<sup>1*</sup> & David J Civitello<sup>1</sup>   

### _<sup>1</sup> Department of Biology, Emory University, 1510 Clifton Road NE, Atlanta, GA, USA, 30322_    

#### *Corresponding author: matthew.malishev@gmail.com    

doi: [DOI](</>)  

## Overview    

Individual-based model (IBM) for schistosome parasite and snail host population interactions based on dynamic energy budget (DEB) and transmission model. Model outlines how resource use by a size-structured host population influences population cycles of hosts and parasites. Host and parasite ecology and energetics are captured using a mechanistic energy budget model and interactions among hosts and parasites are captured by a disease tranmission model. Resource use by hosts and thus parasites is captured by the IBM.       

Files:  

.R    
.html    
.pdf  
.nlogo    
.c  
.o  
.so  
.sh  

******  

## Running the simulation model  

Download the instructions for your operating system  

[*Windows*](https://github.com/darwinanddavis/SchistoIBM/tree/master/windows)  
[*Mac OSX*](https://github.com/darwinanddavis/SchistoIBM/tree/master/mac)  

## :pig: Troubleshooting for running Netlogo from `R` on Mac OSX. See the [Github page](https://github.com/darwinanddavis/rnetlogo_diagnostics) for a detailed breakdown of the troubleshooting steps.  
:one: [Installing compiler toolchain for Mac OSX](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/)    
:two: if rJava error, run the following in terminal (src: https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite and #http://paulklemm.com/blog/2015-02-20-run-rjava-with-rstudio-under-osx-10-dot-10/):    
``` {bash}
sudo ln -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib 
```  

## References  
[Civitello, D. J., Fatima, H. , Johnson, L. R., Nisbet, R. M., Rohr, J. R. and Ben‐Ami, F. (2018), Bioenergetic theory predicts infection dynamics of human schistosomes in intermediate host snails across ecological gradients. Ecol Lett, 21: 692-701](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12937)  

## Maintainer  
Matt Malishev [@darwinanddavis](https://www.researchgate.net/profile/Matt_Malishev)  

