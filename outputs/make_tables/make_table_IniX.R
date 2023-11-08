# ==============================================================================
# Prepare tex table with initial value of X
# ==============================================================================

# Define format of figures:
nb.dec <- 2 # number of decimal numbers
format.nb  <- paste("%.",nb.dec,"f",sep="")
format.nb0 <- paste("%.",0,"f",sep="")
format.nb1 <- paste("%.",1,"f",sep="")
format.nb5 <- paste("%.",5,"f",sep="")

latex.table <- rbind(
  paste("Initial emissions &",
        "$\\mathcal{E}_0$ &",
        "\\eqref{eq:Emit}&",
        make.entry(model$vector.ini$ini_E,format.nb),"&GtCO$_2$ &",
        "DICE2023",
        "\\\\",sep=""),
  paste("Initial industrial emissions &",
        "$\\mathcal{E}_{ind,0}$ &",
        "\\eqref{eq:qt}&",
        make.entry(model$vector.ini$ini_Eind,format.nb),"&GtCO$_2$&",
        "DICE2023",
        "\\\\",sep=""),
  paste("Initial radiative forcings &",
        "$F_0$ &",
        "\\eqref{eq:Fradiat}+\\eqref{eq:Tat}&",
        make.entry(model$vector.ini$ini_F,format.nb),"&W/m$^2$&",
        "CDICE",
        "\\\\",sep=""),
  paste("Initial carbon concent. (atmos.) &",
        "$M_{AT,0}$ &",
        "\\eqref{eq:Mvector}&",
        make.entry(model$vector.ini$ini_Mat,format.nb0),"&GtC&",
        "CDICE",
        "\\\\",sep=""),
  paste("Initial carbon concent. (upper ocean) &",
        "$M_{UP,0}$ &",
        "\\eqref{eq:Mvector}&",
        make.entry(model$vector.ini$ini_Mup,format.nb0),"&GtC&",
        "CDICE",
        "\\\\",sep=""),
  paste("Initial carbon concent. (lower ocean) &",
        "$M_{LO,0}$ &",
        "\\eqref{eq:Mvector}&",
        make.entry(model$vector.ini$ini_Mlo,format.nb0),"&GtC&",
        "CDICE",
        "\\\\",sep=""),
  paste("Initial temp. anomaly (atmos.) &",
        "$T_{AT,0}$ &",
        "\\eqref{eq:Tat}&",
        make.entry(model$vector.ini$ini_Tat,format.nb),"&\\degree C&",
        "CDICE",
        "\\\\",sep=""),
  paste("Initial temp. anomaly (lower ocean) &",
        "$T_{LO,0}$ &",
        "\\eqref{eq:TLO}&",
        make.entry(model$vector.ini$ini_Tlo,format.nb),"&\\degree C&",
        "CDICE",
        "\\\\",sep=""),
  paste("Initial sea level &",
        "$H_{0}$ &",
        "\\eqref{eq:SeaL}&",
        make.entry(model$vector.ini$ini_H,format.nb),"&m&",
        "\\citet{Vermeer_Rahmstorf_2009}",
        "\\\\",sep="")
)

latex.file <- paste("outputs/Tables/table_ini_X.txt",sep="")

write(latex.table, file = latex.file)

