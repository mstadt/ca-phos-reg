caphos_mod <- function(t,y,params) {
    # Variable names
    # FGF23 <- y[1]
    # FGFRbone <- y[2]
    # FGFR <- y[3] # effect of FGF23 on renal PO4 excretion

    dydt <- c()
    with(params, {
        # ECCPhos = ECCPhos0
        # Ctriol = Ctriol0

        # # ddt_FGF23
        # JFGF23 = BFGF23*kFGF
        # StimPhos = (ECCPhos/ECCPhos0)^gamma_phosSTIM
        # StimCtriol = (Ctriol/Ctriol0)^gamma_ctriolSTIM
        # dydt[1] <- JFGF23 * StimPhos * StimCtriol - kFGF23 * FGF23 * FGFRbone

        # # ddt_FGFRbone
        # dydt[2] <- kFGF - kFGF * FGFRbone

        # # ddt_FGFR
        # dydt[3] <- kFGF23 - kFGF23*FGFR


        # list of ODES
        list(c(dydt))
    })
}