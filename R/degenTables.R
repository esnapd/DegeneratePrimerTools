%{$degentables{degen1}} = ( # 1. The Standard Code
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen2}} = ( # 2. The Vertebrate Mitochondrial Code
  #        Code 2          Standard
  # AGA    Ter  *          Arg  R
  # AGG    Ter  *          Arg  R
  # AUA    Met  M          Ile  I
  # UGA    Trp  W          Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen3}} = ( # 3. The Yeast Mitochondrial Code
  #        Code 3          Standard
  # AUA    Met  M          Ile  I
  # CUU    Thr  T          Leu  L
  # CUC    Thr  T          Leu  L
  # CUA    Thr  T          Leu  L
  # CUG    Thr  T          Leu  L
  # UGA    Trp  W          Ter  *
  # CGA    absent          Arg  R
  # CGC    absent          Arg  R
  TTR => [qw(TTA TTG TTR)], # Leu
  TTN => [qw(TTM TTW TTS TTK TTV TTH TTD TTB TTN)], # Leu & Phe
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGG CGR CGY CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen4}} = ( # 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
  #        Code 4         Standard
  # UGA    Trp  W          Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen5}} = ( # 5. The Invertebrate Mitochondrial Code
  #        Code 5          Standard
  # AGA    Ser  S          Arg  R
  # AGG    Ser  S          Arg  R
  # AUA    Met  M          Ile  I
  # UGA    Trp  W          Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen6}} = ( # 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
  #          Code 6       Standard
  # UAA      Gln  Q        Ter  *
  # UAG      Gln  Q        Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen9}} = ( # 9. The Echinoderm and Flatworm Mitochondrial Code
  #        Code 9        Standard
  # AAA      Asn  N        Lys K
  # AGA      Ser  S        Arg R
  # AGG      Ser  S        Arg R
  # UGA      Trp  W        Ter *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen10}} = ( # 10. The Euplotid Nuclear Code
  #          Code 10     Standard
  # UGA      Cys  C        Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen11}} = ( # 11. The Bacterial, Archaeal and Plant Plastid Code
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen12}} = ( # 12. The Alternative Yeast Nuclear Code
  #           Code 12      Standard
  # CUG       Ser          Leu
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu; Phe & Leu; Phe, Leu & Ser3
  CTG => [qw(CTG)], # Ser3
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg; Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degen13}} = ( # 13. The Ascidian Mitochondrial Code
  #        Code 13         Standard
  # AGA    Gly  G          Arg  R
  # AGG    Gly  G          Arg  R
  # AUA    Met  M          Ile  I
  # UGA    Trp  W          Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
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
  RGN => [qw(AGM AGK AGS AGW AGH AGB AGV AGD AGN AGA AGG AGR GGA GGC GGG GGT GGR GGY GGM GGK GGS
             GGW GGH GGB GGV GGD GGN RGC RGA RGT RGG RGR RGY RGM RGK RGS RGW RGH RGB RGV RGD RGN)], # Gly & Ser2
  NNN => [qw(NNN)], # degen
  '---' => [qw(---)], # indel
);

%{$degentables{degen14}} = ( # 14. The Alternative Flatworm Mitochondrial Code
  #          Code 14      Standard
  # AAA      Asn  N       Lys  K
  # AGA      Ser  S       Arg  R
  # AGG      Ser  S       Arg  R
  # UAA      Tyr  Y       Ter  *
  # UGA      Trp  W       Ter  *
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
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
  '---' => [qw(---)], # indel
);

%{$degentables{degenS}} = (
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degenSZ}} = (
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

%{$degentables{degenZ}} = (
  YTN => [qw(TTM TTK TTS TTW TTH TTB TTV TTD TTN TTA TTG TTR CTA CTC CTG CTT CTR CTY CTM CTK CTS
             CTW CTH CTB CTV CTD CTN YTC YTA YTT YTG YTR YTY YTM YTK YTS YTW YTH YTB YTV YTD YTN)], # Leu, Phe & Leu
  TTY => [qw(TTT TTC TTY)], # Phe
  MGN => [qw(CGT CGA CGC CGG CGR CGY CGM CGK CGS CGW CGH CGB CGV CGD CGN AGA AGG AGR AGM AGK AGS
             AGW AGH AGB AGV AGD AGN MGC MGA MGT MGG MGR MGY MGM MGK MGS MGW MGH MGB MGV MGD MGN)], # Arg, Arg & Ser2
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
  '---' => [qw(---)], # indel
);

$tablenumber = 1; # set default table number to 1 (standard genetic code)


# MAIN ######################################################################

&parse_parameters;
print "\n DEGENERACY CODING $version \t(c) Andreas Zwick & April Hussey 2008-2012\n";
print "\************************\n";
print "Genetic code table: '$tablenumber'.\n";
&prep_table;
&get_data;
&degenerate;
&create_output;
print "\nALL DONE!\n\n";

# END MAIN ##################################################################


sub parse_parameters () {
  my ($parameter, $table_name);
  foreach $parameter (@ARGV) {	#parse all parameters
    if (($parameter =~ /^\-t\ *=\ *\d+/i) || ($parameter =~ /^--table\ *=\ *\d+/i)){
      $parameter =~ s/^.+=//;
      $tablenumber = uc($parameter);
      if (! exists $degentables{"degen" . $tablenumber}) {
        &help_text("Don't know the specified table: '$tablenumber'!");
      }
      next;
    } else {
      if ( ! -e $parameter) {
        &help_text("Don't know what to do with the parameter: '$parameter'!");
        next;
      } else {
        $filename = $parameter;
      }
    }
  }
  %degencodons = %{$degentables{"degen" . $tablenumber}};
} #END parse_parameters


sub prep_table () { # check each element in the hash
  foreach $codonkey (keys %degencodons) {
    @codonarray = @{$degencodons{$codonkey}};
    #print "@codonarray, $codonkey\n";
    foreach $codon (@codonarray) {
      if (exists $degencodons1by1{$codon}) {
        $warnings .= "codon ($codon) already exists as a key\n";
        print "codon ($codon) already exists as a key\n";
      }
      $degencodons1by1{$codon}=$codonkey;
    }
  }
} # END prep_table


sub get_data () {
  # Store file contents in an array
  open(INPUT, "$filename") or die "An error occured while attempting to open (R) $filename: $!";
  print "Reading input file '$filename'.\n";
  @filecontents = <INPUT>;
  close(INPUT);

  #replace any spaces and tabs in fasta or flat file sequence identifiers with _
  foreach (@filecontents) {
    if (m/^[\>\#]/i) {
             s/[\ \t]/_/g;
  }
}

# Remove ending from input filename; recognizes filenames with internal dots correctly
$filename =~ s/(.*\.*)(\..*)/$1/;

# Extract path to file, if input file had been specified with a path
if ($filename =~ s/(.*[\/\\])//) {
  $filepath = $1;
}

# Create a directory for the output files, then enter that directory
mkdir("Degen_$filename", 0755);
chdir("Degen_$filename");

# turn matrix in (interleaved) fasta, flat or nexus format into a hash matrix
$concatcontents = join('<', @filecontents);  # join the array into a single string, marking linebreaks with '<'
$concatcontents =~ s/[\#\>]/!/g;             # convert > and # of FLAT and FASTA format into identifier marker !
                       $concatcontents =~ s/!\s+/!/g;               # elimiante any whitespaces preceding the sequence identifiers
                       $concatcontents =~ s/\<\ *\[.*\]\ */!/g;     # turn square brackets at start of line into marker for sequence identifier (e.g., seaview nexus file)
                       $concatcontents =~ s/\ *\[.*\]\ *//g;        # eliminate any text in square brackets (e.g., comments from seaview)
                       $concatcontents =~ s/[\ \t]+/@/g;            # replace any number of spaces or tabs with a single whitespace marker @
                       $concatcontents =~ s/\r?\n|\r//g;            # eliminate any linebreaks
                       $concatcontents =~ s/.*Matrix//i;            # eliminate any NEXUS structures before sequences
                       $concatcontents =~ s/\s*;.*//;               # eliminate any NEXUS structures after sequences
                       $concatcontents =~ s/[\<\!\@][\<\!\@]+/!/g;  # replace multiple subsequent linebreak markers with a single linebreak marker !
                       $concatcontents =~ s/(\!.+?)\</$1!/g;        # enclose sequence identifiers with ! ... !
                       $concatcontents =~ s/\<(.+?)\@/!$1!/g;       # enclose sequence identifiers with ! ... !
                       $concatcontents =~ s/\<//g;                  # remove unnecessary linebreak markers
                       @sequences = split(/[\<\>\@\!]/, $concatcontents);
                       shift (@sequences);
                       for ($i = 0; $i < $#sequences; $i=$i+2) {
                            $hashofsequences{$sequences[$i]} .= $sequences[$i+1];
                            }
