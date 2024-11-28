

model_Traeger <- model_sol

model_Traeger$target_vector["ECumD2"]   <- compute_alternative_damage(2,type="L. Traeger")
model_Traeger$target_vector["ECumD4"]   <- compute_alternative_damage(4,type="L. Traeger")
model_Traeger$target_vector["stdCumD4"] <- .000001
#model_Traeger$target_vector["mu_c0"] <- .08
model_Traeger$target_vector["sigma_c0"] <- .000001

model_Traeger <- solveParam4D(model_Traeger)


#model_Traeger$parameters$a_D <- model_Traeger$parameters$a_D/3
#model_Traeger$parameters$b_D <- model_Traeger$parameters$b_D/3
model_Traeger$parameters$a_N <- 0.00000001
model_Traeger$parameters$b_N <- 0.00000001
model_Traeger$parameters$mu_D <- 0.00000001
model_Traeger$parameters$mu_N <- 0.00000001
model_Traeger$parameters$mu_T <- 0.00000001
model_Traeger$parameters$mu_H <- 0.00000001
#model_Traeger$parameters$sigma_a <- 0.00000001
model_Traeger$parameters$b_sk <- 0

model_Traeger$parameters$gamma <- 1.0001

model_Traeger$parameters$delta <- (1 - .014)^5
#model_Traeger$parameters$delta <- (1 - .005)^5
model_Traeger <- solveParam4c(model_Traeger)

#model_Traeger$parameters$delta_K <- .4
#model_Traeger$parameters$A_bar <- .2

model_sol_Traeger <- model_solve(model_Traeger,indic_mitig = TRUE,
                                 indic_CRRA = FALSE)

# model_sol_Traeger_2 <- model_solve(model_Traeger,indic_mitig = FALSE,
#                                    mu.chosen = 1*pmax(model_sol_Traeger$mu,1),
#                                    indic_CRRA = FALSE)


SCC <- scc.fct(model_sol_Traeger,h=0)
print(SCC)

EV <- EV.fct(model_sol_Traeger)
plot(EV$EX$T_at)
plot(EV$EX$delc)
plot(model_sol_Traeger$mu)

omega_ZCB <- matrix(0,model_sol$n.X)
Price.ZC <- varphi(model_sol_Traeger,
                   omega.varphi=omega_ZCB,
                   H = 30)
yds_CR_EZ <- Price.ZC$r.t/100
plot(yds_CR_EZ,type="l")


