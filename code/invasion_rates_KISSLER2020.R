# =============================================================================
# Copyright
# =============================================================================
# COPYRIGHT 2020 STEPHEN M. KISSLER
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.


# =============================================================================
# Internal functions
# =============================================================================
seir_model_2strains <- function(t, y, parms) {
  dS1S2 <- -S1S2c1(t,y,parms) - S1S2c2(t,y,parms) + 
    R1S2c1(t,y,parms) + S1R2c2(t,y,parms)
  dE1S2 <- -E1S2c1(t,y,parms) - E1S2c2(t,y,parms) + 
    S1S2c1(t,y,parms) + E1R2c2(t,y,parms)
  dS1E2 <- -S1E2c1(t,y,parms) - S1E2c2(t,y,parms) + 
    R1E2c1(t,y,parms) + S1S2c2(t,y,parms)
  dE1E2 <- -E1E2c1(t,y,parms) - E1E2c2(t,y,parms) + 
    S1E2c1(t,y,parms) + E1S2c2(t,y,parms)
  dI1S2 <- -I1S2c1(t,y,parms) - I1S2c2(t,y,parms) + 
    E1S2c1(t,y,parms) + I1R2c2(t,y,parms)
  dS1I2 <- -S1I2c1(t,y,parms) - S1I2c2(t,y,parms) + 
    R1I2c1(t,y,parms) + S1E2c2(t,y,parms)
  dR1S2 <- -R1S2c1(t,y,parms) - R1S2c2(t,y,parms) + 
    I1S2c1(t,y,parms) + R1R2c2(t,y,parms)
  dI1E2 <- -I1E2c1(t,y,parms) - I1E2c2(t,y,parms) + 
    E1E2c1(t,y,parms) + I1S2c2(t,y,parms)
  dE1I2 <- -E1I2c1(t,y,parms) - E1I2c2(t,y,parms) + 
    S1I2c1(t,y,parms) + E1E2c2(t,y,parms)
  dS1R2 <- -S1R2c1(t,y,parms) - S1R2c2(t,y,parms) + 
    R1R2c1(t,y,parms) + S1I2c2(t,y,parms)
  dR1E2 <- -R1E2c1(t,y,parms) - R1E2c2(t,y,parms) + 
    I1E2c1(t,y,parms) + R1S2c2(t,y,parms)
  dI1I2 <- -I1I2c1(t,y,parms) - I1I2c2(t,y,parms) + 
    E1I2c1(t,y,parms) + I1E2c2(t,y,parms)
  dE1R2 <- -E1R2c1(t,y,parms) - E1R2c2(t,y,parms) + 
    S1R2c1(t,y,parms) + E1I2c2(t,y,parms)
  dR1I2 <- -R1I2c1(t,y,parms) - R1I2c2(t,y,parms) + 
    I1I2c1(t,y,parms) + R1E2c2(t,y,parms)
  dI1R2 <- -I1R2c1(t,y,parms) - I1R2c2(t,y,parms) + 
    E1R2c1(t,y,parms) + I1I2c2(t,y,parms)
  dR1R2 <- -R1R2c1(t,y,parms) - R1R2c2(t,y,parms) + 
    I1R2c1(t,y,parms) + R1I2c2(t,y,parms)
  dinc1 <- S1S2c1(t,y,parms) + S1E2c1(t,y,parms) + S1I2c1(t,y,parms) + S1R2c1(t,y,parms)
  dinc2 <- S1S2c2(t,y,parms) + E1S2c2(t,y,parms) + I1S2c2(t,y,parms) + R1S2c2(t,y,parms)
  
  return(list(c(dS1S2,dE1S2,dS1E2,dE1E2,dI1S2,dS1I2,dR1S2,dI1E2,dE1I2,dS1R2,dR1E2,dI1I2,dE1R2,dR1I2,dI1R2,dR1R2,dinc1,dinc2)))
}

beta.val <- function(t, amplitude=1, baseline=1.5, phi.val=-4, gamma.val=1,tstep=1/7){
	gamma.val*(amplitude/2 * cos(2*pi*(t-phi.val)/52*tstep) + (amplitude/2 + baseline))
}



import_strain1 <- function(t,kappa.val,importtime1,importlength){
	ifelse(t>importtime1 & t<=(importtime1+importlength), kappa.val, 0)
}

import_strain2 <- function(t,kappa.val,importtime2,importlength){
	ifelse(t>importtime2 & t<=(importtime2+importlength), kappa.val, 0)
}

# =============================================================================
# State update functions 
# =============================================================================

S1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1S2 + import_strain1(t,kappa.val,importtime1,importlength)*S1S2
  beta.val1*(I1S2+I1E2+I1I2+I1R2)*S1S2 + import_strain1(t,kappa.val,importtime1,importlength)*S1S2
	})}
S1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*S1S2 + import_strain2(t,kappa.val,importtime2,importlength)*S1S2
  beta.val2*(S1I2+E1I2+I1I2+R1I2)*S1S2 + import_strain2(t,kappa.val,importtime2,importlength)*S1S2
	})}
E1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1S2
	})}
E1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*E1S2 + (1-chi12.val)*import_strain2(t,kappa.val,importtime2,importlength)*E1S2
  (1-chi12.val)*beta.val2*(S1I2+E1I2+I1I2+R1I2)*E1S2 + (1-chi12.val)*import_strain2(t,kappa.val,importtime2,importlength)*E1S2
	})}
S1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1E2 + (1-chi21.val)*import_strain1(t,kappa.val,importtime1,importlength)*S1E2
  (1-chi21.val)*beta.val1*(I1S2+I1E2+I1I2+I1R2)*S1E2 + (1-chi21.val)*import_strain1(t,kappa.val,importtime1,importlength)*S1E2
  
	})}
S1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*S1E2
	})}
E1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1E2 
	})}
E1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1E2 
	})}
I1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1S2
	})}
I1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*I1S2 + (1-chi12.val)*import_strain2(t,kappa.val,importtime2,importlength)*I1S2
  (1-chi12.val)*beta.val2*(S1I2+E1I2+I1I2+R1I2)*I1S2 + (1-chi12.val)*import_strain2(t,kappa.val,importtime2,importlength)*I1S2
  
	})}
S1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1I2 + (1-chi21.val)*import_strain1(t,kappa.val,importtime1,importlength)*S1I2
  (1-chi21.val)*beta.val1*(I1S2+I1E2+I1I2+I1R2)*S1I2 + (1-chi21.val)*import_strain1(t,kappa.val,importtime1,importlength)*S1I2
  
	})}
S1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*S1I2
	})}
R1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1S2
	})}
R1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*R1S2 + (1-chi12.val)*import_strain2(t,kappa.val,importtime2,importlength)*R1S2
  (1-chi12.val)*beta.val2*(S1I2+E1I2+I1I2+R1I2)*R1S2 + (1-chi12.val)*import_strain2(t,kappa.val,importtime2,importlength)*R1S2
  
	})}
I1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1E2
	})}
I1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*I1E2
	})}
E1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1I2
	})}
E1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*E1I2
	})}
S1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	#(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1R2 + (1-chi21.val)*import_strain1(t,kappa.val,importtime1,importlength)*S1R2
  (1-chi21.val)*beta.val1*(I1S2+I1E2+I1I2+I1R2)*S1R2 + (1-chi21.val)*import_strain1(t,kappa.val,importtime1,importlength)*S1R2
  
	})}
S1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*S1R2
	})}
R1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1E2
	})}
R1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*R1E2
	})}
I1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1I2 
	})}
I1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1I2
	})}
E1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1R2 
	})}
E1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*E1R2
	})}
R1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1I2
	})}
R1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*R1I2
	})}
I1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1R2
	})}
I1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*I1R2
	})}
R1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1R2
	})}
R1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*R1R2
	})}
