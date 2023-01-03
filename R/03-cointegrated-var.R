
source("./R/01-data.R")



source("https://raw.githubusercontent.com/mmoessler/juselius-2006/main/R/ca_jo_jus06_fun.R")



# with trend in cir, with shift dummy and with Ds831 in cir
ca.jo.res.01 <- ca_jo_jus06_fun(x = data[,c("Lm3rC","Lyr","Dpy","Rm","Rb")], type = c("trace"), ecdet = c("trend"), K = 2, spec = c("transitory"), data = data)



# without trend in cir, with shift dummy and with Ds831 in cir
ca.jo.res.02 <- ca_jo_jus06_fun(x = data[,c("Lm3rC","Lyr","Dpy","Rm","Rb")], type = c("trace"), ecdet = c("none"), K = 2, spec = c("transitory"), data = data)



round(Re(ca.jo.res.01$lambda), 2) 



round(Re(ca.jo.res.01$Vorg), 2)



V.01 <- cbind(cbind(ca.jo.res.01$Vorg[,1]/ca.jo.res.01$Vorg[3,1]),
              cbind(ca.jo.res.01$Vorg[,2]/ca.jo.res.01$Vorg[1,2]),
              cbind(ca.jo.res.01$Vorg[,3]/ca.jo.res.01$Vorg[4,3]),
              cbind(ca.jo.res.01$Vorg[,4]/ca.jo.res.01$Vorg[2,4]),
              cbind(ca.jo.res.01$Vorg[,5]/ca.jo.res.01$Vorg[5,5]),
              cbind(ca.jo.res.01$Vorg[,6]/ca.jo.res.01$Vorg[1,6]),
              cbind(ca.jo.res.01$Vorg[,7]/ca.jo.res.01$Vorg[1,7]))

round(V.01, 2)



W.01 <- ca.jo.res.01$S0K %*% V.01 %*% solve(t(V.01) %*% ca.jo.res.01$SKK %*% V.01)

round(W.01, 2)



round(ca.jo.res.01$PI, 2)

