#!/usr/bin/env python

# Create Degeneracy Tables Using The Degen.pl Codon Substitution Tables
# http://www.phylotools.com/ptdegendocumentation.htm
# Data from the dDegen_v1_4.pl was copied, pasted, and lightly organized to
# facilitate parsing

import operator
from collections import OrderedDict

# 1. The Standard Code
degen1 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], #  # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 2. The Vertebrate Mitochondrial Code
#        Code 2          Standard
# AGA    Ter  *          Arg  R
# AGG    Ter  *          Arg  R
# AUA    Met  M          Ile  I
# UGA    Trp  W          Ter  *
degen2 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
CGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN)], # Arg
AGY => [qw(AGC AGT AGY)], # Ser2
ATY => [qw(ATT ATC ATY)], # Ile
ATN => [qw(ATM ATW ATS ATK ATV ATH ATD ATB ATN)], # Ile & Met
ATR => [qw(ATG ATA ATR)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGN => [qw(TGM TGW TGS TGK TGV TGH TGD TGB TGN)], # Cys & Trp
TGR => [qw(TGG TGA TGR)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 3. The Yeast Mitochondrial Code
#        Code 3          Standard
# AUA    Met  M          Ile  I
# CUU    Thr  T          Leu  L
# CUC    Thr  T          Leu  L
# CUA    Thr  T          Leu  L
# CUG    Thr  T          Leu  L
# UGA    Trp  W          Ter  *
# CGA    absent          Arg  R
# CGC    absent          Arg  R
degen3 = """TTR => [qw(TTA TTG TTR)], # Leu
TTN => [qw(TTM TTW TTS TTK TTV TTH TTD TTB TTN)], # Leu & Phe
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGG CGR CGY CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATY => [qw(ATT ATC ATY)], # Ile
ATN => [qw(ATM ATW ATS ATK ATV ATH ATD ATB ATN)], # Ile & Met
ATR => [qw(ATG ATA ATR)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr1
CTN => [qw(CTT CTA CTC CTG CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN)], # Thr2
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGN => [qw(TGM TGW TGS TGK TGV TGH TGD TGB TGN)], # Cys & Trp
TGR => [qw(TGG TGA TGR)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
#        Code 4         Standard
# UGA    Trp  W          Ter  *
degen4  = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGM TGK TGS TGW TGH TGB TGV TGD TGN)], # Cys & Trp
TGR => [qw(TGA TGG TGR)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 5. The Invertebrate Mitochondrial Code
#        Code 5          Standard
# AGA    Ser  S          Arg  R
# AGG    Ser  S          Arg  R
# AUA    Met  M          Ile  I
# UGA    Trp  W          Ter  *
degen5 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
CGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN)], # Arg
AGN => [qw(AGT AGA AGC AGG AGR AGY AGM AGK AGS AGW AGH AGB AGV AGD AGN)], # Ser2
ATY => [qw(ATT ATC ATY)], # Ile
ATN => [qw(ATM ATW ATS ATK ATV ATH ATD ATB ATN)], # Ile & Met
ATR => [qw(ATG ATA ATR)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGN => [qw(TGM TGW TGS TGK TGV TGH TGD TGB TGN)], # Cys & Trp
TGR => [qw(TGG TGA TGR)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
#          Code 6       Standard
# UAA      Gln  Q        Ter  *
# UAG      Gln  Q        Ter  *
degen6 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
YAR => [qw(TAA TAG TAR CAA CAG CAR YAA YAG YAR)], # Gln
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
YAN => [qw(TAM TAK TAS TAW TAH TAB TAV TAD TAN CAM CAK CAS CAW CAH CAB CAV CAD CAN YAM YAK YAS YAW YAH YAB YAV YAD YAN)], # Gln, Tyr & His
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 9. The Echinoderm and Flatworm Mitochondrial Code
#        Code 9        Standard
# AAA      Asn  N        Lys K
# AGA      Ser  S        Arg R
# AGG      Ser  S        Arg R
# UGA      Trp  W        Ter *
degen9 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
CGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN)], # Arg
AGN => [qw(AGT AGA AGC AGG AGR AGY AGM AGK AGS AGW AGH AGB AGV AGD AGN)], # Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAH => [qw(AAA AAT AAC AAY AAM AAW AAH)], # Asn
AAN => [qw(AAR AAK AAS AAB AAV AAD AAN)], # Asn & Lys
AAG => [qw(AAG)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGN => [qw(TGM TGW TGS TGK TGV TGH TGD TGB TGN)], # Cys & Trp
TGR => [qw(TGG TGA TGR)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 10. The Euplotid Nuclear Code
#          Code 10     Standard
# UGA      Cys  C        Ter  *
degen10 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGH => [qw(TGT TGC TGA TGH TGY TGM TGW)], # Cys
TGN => [qw(TGR TGK TGS TGB TGV TGD TGN)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 11. The Bacterial, Archaeal and Plant Plastid Code
degen11 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 12. The Alternative Yeast Nuclear Code
#           Code 12      Standard
# CUG       Ser          Leu
degen12 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu; Phe & Leu; Phe, Leu & Ser3
CTG => [qw(CTG)], # Ser3
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg; Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
AGY => [qw(AGT AGC AGY)], # Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 13. The Ascidian Mitochondrial Code
#        Code 13         Standard
# AGA    Gly  G          Arg  R
# AGG    Gly  G          Arg  R
# AUA    Met  M          Ile  I
# UGA    Trp  W          Ter  *
degen13 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
CGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN)], # Arg
AGY => [qw(AGC AGT AGY)], # Ser2
ATY => [qw(ATT ATC ATY)], # Ile
ATN => [qw(ATM ATW ATS ATK ATV ATH ATD ATB ATN)], # Ile & Met
ATR => [qw(ATG ATA ATR)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGN => [qw(TGM TGW TGS TGK TGV TGH TGD TGB TGN)], # Cys & Trp
TGR => [qw(TGG TGA TGR)], # Trp
RGN => [qw(AGM AGK AGS AGW AGH AGB AGV AGD AGN AGA AGG AGR GGA GGC GGG GGT GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN RGC RGA RGT RGG RGR RGY RGM RGK RGS RGW RGH RGB RGV RGD RGN)], # Gly & Ser2
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

# 14. The Alternative Flatworm Mitochondrial Code
#          Code 14      Standard
# AAA      Asn  N       Lys  K
# AGA      Ser  S       Arg  R
# AGG      Ser  S       Arg  R
# UAA      Tyr  Y       Ter  *
# UGA      Trp  W       Ter  *
degen14 = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
CGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN)], # Arg
AGN => [qw(AGT AGA AGC AGG AGR AGY AGM AGK AGS AGW AGH AGB AGV AGD AGN)], # Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN)], # Ser1
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAH => [qw(TAT TAC TAA TAH TAY TAM TAW)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAH => [qw(AAA AAT AAC AAY AAM AAW AAH)], # Asn
AAN => [qw(AAR AAK AAS AAB AAV AAD AAN)], # Asn & Lys
AAG => [qw(AAG)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGN => [qw(TGM TGW TGS TGK TGV TGH TGD TGB TGN)], # Cys & Trp
TGR => [qw(TGG TGA TGR)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

degenS = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
AGY => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN AGT AGC AGY)], # Ser1 and Ser2 to Ser2
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""

degenSZ = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN AGT AGC AGY)], # degen
--- => [qw(---)], # indel
"""

degenZ = """YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
TTY => [qw(TTT TTC TTY)], # Phe
MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
ATH => [qw(ATT ATC ATA ATH ATY ATM ATW)], # Ile
ATN => [qw(ATR ATK ATS ATB ATV ATD ATN)], # Ile & Met
ATG => [qw(ATG)], # Met
GTN => [qw(GTT GTA GTC GTG GTR GTY GTM GTK GTS GTW GTH GTB GTV GTD GTN)], # Val
TCN => [qw(TCT TCA TCC TCG TCR TCY TCM TCK TCS TCW TCH TCB TCV TCD TCN AGT AGC AGY)], # Ser1 and Ser2 to Ser1
CCN => [qw(CCT CCA CCC CCG CCR CCY CCM CCK CCS CCW CCH CCB CCV CCD CCN)], # Pro
ACN => [qw(ACT ACA ACC ACG ACR ACY ACM ACK ACS ACW ACH ACB ACV ACD ACN)], # Thr
GCN => [qw(GCT GCA GCC GCG GCR GCY GCM GCK GCS GCW GCH GCB GCV GCD GCN)], # Ala
TAY => [qw(TAT TAC TAY)], # Tyr
CAY => [qw(CAT CAC CAY)], # His
CAN => [qw(CAM CAK CAS CAW CAH CAB CAV CAD CAN)], # His & Gln
CAR => [qw(CAA CAG CAR)], # Gln
AAY => [qw(AAT AAC AAY)], # Asn
AAN => [qw(AAM AAK AAS AAW AAH AAB AAV AAD AAN)], # Asn & Lys
AAR => [qw(AAA AAG AAR)], # Lys
GAY => [qw(GAT GAC GAY)], # Asp
GAN => [qw(GAM GAK GAS GAW GAH GAB GAV GAD GAN)], # Asp & Glu
GAR => [qw(GAA GAG GAR)], # Glu
TGY => [qw(TGT TGC TGY)], # Cys
TGB => [qw(TGK TGS TGB)], # Cys & Trp
TGG => [qw(TGG)], # Trp
GGN => [qw(GGT GGA GGC GGG GGR GGY GGM GGK GGS GGW GGH GGB GGV GGD GGN)], # Gly
NNN => [qw(NNN)], # degen
--- => [qw(---)], # indel
"""



degentables = {"degen1":  degen1,  "degen2":  degen2,  "degen3":  degen3,  "degen4":  degen4, "degen5":  degen5,
               "degen6":  degen6,  "degen9":  degen9, "degen10": degen10,
               "degen11": degen11, "degen12": degen12, "degen13": degen13, "degen14": degen14,
               "degenS":  degenS,  "degenSZ": degenSZ, "degenZ":  degenZ}

def make_orderedDict(s):
    od = OrderedDict()
    s1 = s.replace("=> [qw(", ":")
    s1 = s1.replace(")],", ":")
    s1 = s1.replace(") ],", ":")

    for line in s1.splitlines():
        target, values, AAs = line.split(":")
        #print(target, values)
        for value in values.split():
            #print(value.strip())
            od[value.strip()] = target.strip()
    return od

def convert_to_R(d, name, prnt=False):
    sorted_d = sorted(d.items(), key=operator.itemgetter(0))
    table_vals =  []
    nameline = "{} <- list(".format(name)
    table_vals.append(nameline)
    if prnt: print(nameline)

    for idx, (k, v) in enumerate(sorted_d):
        if idx == len(sorted_d)-1:
            lookup = '"{}"="{}")'.format(k,v)
            #print(idx, lookup, k, v)
        else:
            lookup = '"{}"="{}",'.format(k,v)
        table_vals.append(lookup)
        if prnt: print(lookup)

    return("\n".join(table_vals))

with open("degeneracy_tables.R","w") as f:
    for k,v in degentables.items():
        table = make_orderedDict(v)
        R_table = convert_to_R(table, k)
        f.write(R_table)
        f.write("\n")


