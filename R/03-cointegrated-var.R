
# Load data
source("./R/01-data.R")

source("https://raw.githubusercontent.com/mmoessler/juselius-2006/main/R/ca_jo_jus06_fun.R")

# with trend in cir, with shift dummy and with Ds831 in cir
ca.jo.res.01 <- ca_jo_jus06_fun(x = data[,c("Lm3rC","Lyr","Dpy","Rm","Rb")], type = c("trace"), ecdet = c("trend"), K = 2, spec = c("transitory"), data = data)

# without trend in cir, with shift dummy and with Ds831 in cir
ca.jo.res.02 <- ca_jo_jus06_fun(x = data[,c("Lm3rC","Lyr","Dpy","Rm","Rb")], type = c("trace"), ecdet = c("none"), K = 2, spec = c("transitory"), data = data)
