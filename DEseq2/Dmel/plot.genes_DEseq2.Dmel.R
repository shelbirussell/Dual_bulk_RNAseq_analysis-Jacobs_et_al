###############################
###### PLOTTING FUNCTION ######

plot_gene <- function(dds1, ID, gene) {

pdf1 <- file.path(paste0("cell_KallistoDESeq2_", ID, "_", gene, "_counts.pdf"))
pdflog2 <- file.path(paste0("cell_KallistoDESeq2_", ID, "_", gene, "_log2.counts.pdf"))

obj <- plotCounts(dds1, gene=ID, intgroup="group", returnData=TRUE)
pdf(pdf1)
print(ggplot(data=obj, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal())
dev.off()
pdf(pdflog2)
print(ggplot(data=obj, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal(panel.grid.minor = element_blank()))
dev.off()
}


#######################
## Module 11 - JW18DOX up / JW18wMel down
####
####
plot_gene(dds, "Dmel_CG10489", "DNApol-epsilon58")
plot_gene(dds, "Dmel_CG10890", "Mus201")
plot_gene(dds, "Dmel_CG1119", "Gnf1")
plot_gene(dds, "Dmel_CG8153", "mus210")
plot_gene(dds, "Dmel_CG6768", "DNApol-epsilon255")

###########
## Module 3 - JW18DOX up
### spliceosome
plot_gene(dds, "Dmel_CG10203", "x16")
plot_gene(dds, "Dmel_CG18591", "SmE")
plot_gene(dds, "Dmel_CG4980", "BCAS2") ## maybe not DE
plot_gene(dds, "Dmel_CG6841", "Prp6")
plot_gene(dds, "Dmel_CG4886", "cyp33")
plot_gene(dds, "Dmel_CG6905", "Cdc5")
plot_gene(dds, "Dmel_CG3605", "Sf3b2") ## not DE
plot_gene(dds, "Dmel_CG5931", "l(3)72Ab") ## aka. Brr2
plot_gene(dds, "Dmel_CG5454", "snRNP-U1-C")
plot_gene(dds, "Dmel_CG9998", "U2af50")
plot_gene(dds, "Dmel_CG8749", "snRNP70K")
plot_gene(dds, "Dmel_CG4849", "GTPase.U5snRNAbinding.tef")

############ Module 3 MMR
plot_gene(dds, "Dmel_CG10387", "tos")
plot_gene(dds, "Dmel_CG7003", "Msh6")
plot_gene(dds, "Dmel_CG11482", "Mlh1")
plot_gene(dds, "Dmel_CG6258", "RfC38")
plot_gene(dds, "Dmel_CG9633", "RpA-70")
plot_gene(dds, "Dmel_CG8142", "Part.of.Elg1.RFC-like.complex.glial")
plot_gene(dds, "Dmel_CG5949", "POLD")
plot_gene(dds, "Dmel_CG9193", "PCNA")
plot_gene(dds, "Dmel_CG5602", "DNA−lig1")

## NER
plot_gene(dds, "Dmel_CG10215", "Ercc1")
plot_gene(dds, "Dmel_CG8019", "hay")

## BER
plot_gene(dds, "Dmel_CG13399", "Chrac-14") ## Chromatin accessibility complex 14kD protein (Chrac-14) encodes a subunit of the chromatin accessibility complex, which is an evolutionary conserved nucleosome sliding enzyme. It is also important for DNA damage response
plot_gene(dds, "Dmel_CG1981", "Thd1") ## Enables G/U mismatch-specific uracil-DNA glycosylase activity and double-stranded DNA binding activity
plot_gene(dds, "Dmel_CG12223", "Dsp1") ## Binds preferentially single-stranded DNA and unwinds double-stranded DNA

## DNA replication
plot_gene(dds, "Dmel_CG5553", "S(faf)240")
plot_gene(dds, "Dmel_CG4039", "Mcm6")
plot_gene(dds, "Dmel_CG7108", "Pri50") #together with the subunit encoded by Prim2, is responsible for primase activity (synthesis of short RNA primers). This activity is required during DNA replication ... Okazaki fragment synthesis.
plot_gene(dds, "Dmel_CG4978", "Mcm7")
plot_gene(dds, "Dmel_CG4978", "dpa")
plot_gene(dds, "Dmel_CG5923", "PolA2")
plot_gene(dds, "Dmel_CG5949", "POLA")
plot_gene(dds, "Dmel_CG7538", "Mcm2")
plot_gene(dds, "Dmel_CG4206", "Mcm3")
plot_gene(dds, "Dmel_CG11164", "part.of.ribonuclease.H2.complex")
plot_gene(dds, "Dmel_CG13690", "RNaseHII.subunit.degrades.RNA.of.RNA:DNA")
plot_gene(dds, "Dmel_CG4082", "Mcm5")

###################################
### Module 2
##### inositol phosphate metabolism
plot_gene(dds, "Dmel_CG17028", "CG17028")
plot_gene(dds, "Dmel_CG13688", "Ipk2")
plot_gene(dds, "Dmel_CG3028", "Ipp")
plot_gene(dds, "Dmel_CG17840", "FIG4")
plot_gene(dds, "Dmel_CG2171", "Tpi")
plot_gene(dds, "Dmel_CG9115", "mtm") #  together with the products of Sbf and Rab21, regulates macrophage protrusion formation
plot_gene(dds, "Dmel_CG9784", "phosphoinositide.5-phosphatase") ##  Predicted to be active in neuron projection https://flybase.org/reports/FBgn0030761
plot_gene(dds, "Dmel_CG10260", "PI4KIIIalpha")
plot_gene(dds, "Dmel_CG17027", "inositol.monophosphate1-phosphatase ")
plot_gene(dds, "Dmel_CG42283", "5Ptasel")
plot_gene(dds, "Dmel_CG17029", "CG17029")
plot_gene(dds, "Dmel_CG17471", "PIP4K") # implicated in the regulation of mTOR signalling and control of cell size
plot_gene(dds, "Dmel_CG9245", "Pis")
plot_gene(dds, "Dmel_CG42271", "dmINPP4A")
plot_gene(dds, "Dmel_CG30295", "Ipk1")
plot_gene(dds, "Dmel_CG4317", "Mipp2")
plot_gene(dds, "Dmel_CG10426", "CG10426") ##By controlling the phosphoinositide composition of the cilia membrane of auditory receptor neurons, regulates the cilia localization of ktub and consequently the transient receptor potential channels iav and nompC

####### anemia pathway
plot_gene(dds, "Dmel_CG3697", "mei-9")
plot_gene(dds, "Dmel_", "TOP3beta")
plot_gene(dds, "Dmel_", "FANCI")
plot_gene(dds, "Dmel_", "mus205")
plot_gene(dds, "Dmel_", "DNApol-eta")
plot_gene(dds, "Dmel_", "mus81")
plot_gene(dds, "Dmel_", "spn-A")
plot_gene(dds, "Dmel_", "mms4")
plot_gene(dds, "Dmel_", "Pms2")
plot_gene(dds, "Dmel_", "CG7602")
plot_gene(dds, "Dmel_", "Fancd2")

####### apoptosis - fly
plot_gene(dds, "Dmel_", "mei-41")
plot_gene(dds, "Dmel_", "sav")
plot_gene(dds, "Dmel_", "wts")
plot_gene(dds, "Dmel_", "wgn")
plot_gene(dds, "Dmel_", "hep")
plot_gene(dds, "Dmel_", "Marf")
plot_gene(dds, "Dmel_", "Tspo")
plot_gene(dds, "Dmel_", "Traf6")
plot_gene(dds, "Dmel_", "Drp1")
plot_gene(dds, "Dmel_", "HtrA2")
plot_gene(dds, "Dmel_", "puc")
plot_gene(dds, "Dmel_", "Dcp-1")
plot_gene(dds, "Dmel_", "bic")
plot_gene(dds, "Dmel_", "Bruce")
plot_gene(dds, "Dmel_", "Parp")
plot_gene(dds, "Dmel_", "Jra")
plot_gene(dds, "Dmel_", "rl")
plot_gene(dds, "Dmel_", "crc")
plot_gene(dds, "Dmel_", "tefu")
plot_gene(dds, "Dmel_", "LamC")
plot_gene(dds, "Dmel_", "Pk92B")
plot_gene(dds, "Dmel_", "PEK")
plot_gene(dds, "Dmel_", "Eip74EF")



###################################
### Module 1
## ME1 ribosome biogenesis
plot_gene(dds, "Dmel_", "Ns3")
plot_gene(dds, "Dmel_", "CG3071")
plot_gene(dds, "Dmel_", "CG13097")
plot_gene(dds, "Dmel_", "Rpp30")
plot_gene(dds, "Dmel_", "CG30349")
plot_gene(dds, "Dmel_", "CG7637")
plot_gene(dds, "Dmel_", "CG8368")
plot_gene(dds, "Dmel_", "Nmd3")
plot_gene(dds, "Dmel_", "CG2972")
plot_gene(dds, "Dmel_", "CG3527")
plot_gene(dds, "Dmel_", "Ns1")
plot_gene(dds, "Dmel_", "hoip")
plot_gene(dds, "Dmel_", "CG12301")
plot_gene(dds, "Dmel_", "l(3)72Dn")
plot_gene(dds, "Dmel_", "Nop56")
plot_gene(dds, "Dmel_", "CG4038")
plot_gene(dds, "Dmel_", "CG10314")
plot_gene(dds, "Dmel_", "RIOK2")
plot_gene(dds, "Dmel_", "Fib")
plot_gene(dds, "Dmel_", "l(1)G0020")
plot_gene(dds, "Dmel_", "CG7246")
plot_gene(dds, "Dmel_", "Bka")
plot_gene(dds, "Dmel_", "Ak6")
plot_gene(dds, "Dmel_", "CG11920")
plot_gene(dds, "Dmel_", "CG4806")
plot_gene(dds, "Dmel_", "Nop60B")
plot_gene(dds, "Dmel_", "RIOK1")
plot_gene(dds, "Dmel_", "Rtc1")
plot_gene(dds, "Dmel_", "l(1)G0045")
plot_gene(dds, "Dmel_", "CG12325")
plot_gene(dds, "Dmel_", "CG13185")
plot_gene(dds, "Dmel_", "NHP2")

## ME1 endocytosis
plot_gene(dds, "Dmel_", "step")
plot_gene(dds, "Dmel_", "shi")
plot_gene(dds, "Dmel_", "Vps26")
plot_gene(dds, "Dmel_", "Chc")
plot_gene(dds, "Dmel_", "Rbpn-5")
plot_gene(dds, "Dmel_", "Vps35")
plot_gene(dds, "Dmel_", "Past1")
plot_gene(dds, "Dmel_", "CG2224")
plot_gene(dds, "Dmel_", "Asap")
plot_gene(dds, "Dmel_", "Arpc3B")
plot_gene(dds, "Dmel_", "CHMP2B")
plot_gene(dds, "Dmel_", "Arf79F")
plot_gene(dds, "Dmel_", "Arpc3A")
plot_gene(dds, "Dmel_", "Su(dx)")
plot_gene(dds, "Dmel_", "Arpc2")
plot_gene(dds, "Dmel_", "Rip11")
plot_gene(dds, "Dmel_", "AP-2alpha")
plot_gene(dds, "Dmel_", "mod(r)")
plot_gene(dds, "Dmel_", "Arf51F")
plot_gene(dds, "Dmel_", "Src64B")
plot_gene(dds, "Dmel_", "Snx1")
plot_gene(dds, "Dmel_", "Vrp1")
plot_gene(dds, "Dmel_", "Rab11")
plot_gene(dds, "Dmel_", "Rab7")
plot_gene(dds, "Dmel_", "Vps29")
plot_gene(dds, "Dmel_", "Vps45")
plot_gene(dds, "Dmel_", "sktl")
plot_gene(dds, "Dmel_", "Sara")
plot_gene(dds, "Dmel_", "Khc") ## I already plotted this elsewhere...
plot_gene(dds, "Dmel_", "krz")
plot_gene(dds, "Dmel_", "AP-2mu")
plot_gene(dds, "Dmel_", "Rab10")
plot_gene(dds, "Dmel_", "SWIP")
plot_gene(dds, "Dmel_", "Rab35")
plot_gene(dds, "Dmel_", "Snx3")
plot_gene(dds, "Dmel_", "Hsp70Ab")
plot_gene(dds, "Dmel_", "TSG101")
plot_gene(dds, "Dmel_", "Vps20")


## ME1 ribosome
plot_gene(dds, "Dmel_", "RpL37A")
plot_gene(dds, "Dmel_", "mRpL1")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")

plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "RpS15")
plot_gene(dds, "Dmel_", ""mRpS18C)
plot_gene(dds, "Dmel_", "mRpS6")
plot_gene(dds, "Dmel_", "mRpL13")
plot_gene(dds, "Dmel_", "RpL24-like")
plot_gene(dds, "Dmel_", "mRpS7")
plot_gene(dds, "Dmel_", "RpL39")
plot_gene(dds, "Dmel_", "RpS19b")
plot_gene(dds, "Dmel_", "mRpL20")
plot_gene(dds, "Dmel_", "RpL7")
plot_gene(dds, "Dmel_", "RpS16")
plot_gene(dds, "Dmel_", "mRpL33")
plot_gene(dds, "Dmel_", "mRpS17")
plot_gene(dds, "Dmel_", "RpS28b")
plot_gene(dds, "Dmel_", "RpS12")
plot_gene(dds, "Dmel_", "mRpL16")



plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")
plot_gene(dds, "Dmel_", "")


## ME 10 infection - mito + ecoli homologous ribosomal proteins
plot_gene(dds, "Dmel_CG12373", "mRpL18")
plot_gene(dds, "Dmel_CG7490", "RpLP0")
plot_gene(dds, "Dmel_", "RpL7A")
plot_gene(dds, "Dmel_", "RpS6")
plot_gene(dds, "Dmel_", "RpS15Aa")
plot_gene(dds, "Dmel_", "RpL24")
plot_gene(dds, "Dmel_", "RpS15Ab")
plot_gene(dds, "Dmel_", "RpL3")
plot_gene(dds, "Dmel_", "mRpL23")
