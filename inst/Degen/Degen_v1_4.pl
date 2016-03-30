#!/usr/bin/perl

# Authors: Andreas Zwick & April Hussey
# Created: 08-26-08
# Last Revised: 12-19-12
# Version: 1.4
#
# Copyright (C) 2008-2012  Andreas Zwick & April Hussey
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################
#
# First published as online Appendix S5 in:
# Regier, J.C., Shultz, J.W., Zwick, A., Hussey, A., Ball, B., Wetzer, R.,
# Martin, J.W. & Cunningham, C.W. (2010), Arthropod relationships revealed
# by phylogenomic analysis of nuclear protein-coding sequences.
# Nature 463: 1079-1083.
#
# Alternative encodings published in:
# Zwick A., Regier J.C. & Zwickl, D.J. (2012) Resolving Discrepancy between
# Nucleotides and Amino Acids in Deep-Level Arthropod Phylogenomics:
# Differentiating Serine Codons in 21-Amino-Acid Models.
# PLoS ONE 7(11): e47450. doi:10.1371/journal.pone.0047450
#
###########################################################################
#
# Requirements:
#   - functional PERL installation
#     (tested with v5.14.2 built for x86_64-linux-gnu-thread-multi)
#     (tested with v5.10.0 built for MSWin32-x86-multi-thread)
#     (tested with v5.8.6 built for darwin-thread-multi-level, OS X 10.4.11)
#   - input file in non-interleaved NEXUS, FASTA or FLAT format
#
# Optional:
#   - selection of genetic code (standard code is default)
#
# To execute this script, type on the command line:
#   perl degen_v1_4.pl infile [-t= OR --table=##]
#
# The original input file will not be overwritten.
# Output directory, Degen_<basefilename>, is created in the directory from
# which the script is run and contains 4 output files:
#   1) Degen_<basefilename>.nex - the final dataset of the transformed
#                                 nucleotide sequences in NEXUS format
#   2) Degen_<basefilename>.fasta - the final dataset of the transformed
#                                   nucleotide sequences in FASTA format
#   3) HashRegEx_<basefilename>.txt - the complete listing of all nucleotide
#                                      positions transformed by the script
#   4) Warnings_<basefilename>.txt - the complete listing of all nucleotide
#                                     positions not transformed by the script
#                                     and warnings about unexpected sequence
#                                     lengths or characters.
#   Positions noted in the files are of the first base of the codon mentioned.
#
# For amino acids that are encoded for by disjunct codon groups, these codon
# groups can either be differentiated or merged; the codons that are affected
# differ between organisms with the genetic code. For example, in the standard
# genetic code serine is encoded by two disjunct codon groups, YTN and MGN:
#   A) Fully degenerate each of the two codon groups separately [Degen1]:
#      - Triplets that encode Leu or Leu & Phe are converted to YTN,
#      - triplets that encode Arg or Arg & Ser2 are converted to MGN.
#   B) Degenerate both of the two codon groups to just either of them
#      [DegenS and DegenZ]:
#      - Triplets that encode Ser1 and Ser2 are converted to AGY [DegenS]
#        or TCN [DegenZ]
#   C) Eliminate both of the two codon groups [DegenSZ]
#      - Triplets that encode Ser1 or Ser2 are converted to NNN
# In addition, the following general rules apply:
#      - Triplets that encode His & Gln, Asn & Lys, or Asp & Glu are converted
#        to CAN, AAN, or GAN, respectively
#      - Other triplets that encode nonsynonymous polymorphisms at nt1 and/or
#        nt2 are converted to NNN
#
###########################################################################
#
# CHANGE HISTORY:
# Version 1.3 [19 NOV 2009]:
# - added non-interleaved NEXUS format as input and output
# - recognizes internal dots of filenames correctly, e.g., ART.80tax.68gn.nex
#   is no longer truncated to ART for the output
# - permits the specification of a path to the input file
# - reduced terminal output to increase speed and show only potentially
#   problematic codon changes; to see all changes carried out, check the
#   HashRegEx_... text file
# - file handling errors somewhat more transparent
# - changed suffix of output files to lower case
#
# Version 1.4 [5 JAN 2013]:
# - added support for genetic codes 2-6 and 9-14 of NCBI and for alternative
#   encodings of the standard genetic code (DegenS, DegenSZ and DegenZ)
# - added support for interleaved files as generated by PAUP*, SeaView for
#   Linux and potentially other software
# - fixed miscount of one in the display / record of the error position
#
###########################################################################

use strict;

my $version = "v1.4";
my (%degentables, %degencodons, %degencodons1by1, %hashofsequences, %degensequences);
my (@codonarray, @filecontents, @sequences);
my ($tablenumber, $codonkey, $codon, $warnings, $filename, $filepath, $concatcontents, $i, $seqlength, $hashregexwarnings, $identifier, $position, $codon1);

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
} # END get_data


sub degenerate () {
  # Find highlighted codons in each sequence
  foreach $identifier (keys %hashofsequences) {
    # Find the length of the sequences (they should all be the same length)
    $seqlength = length($hashofsequences{$identifier});
    $warnings .= "Sequence: $identifier, $seqlength bp\n";
    $hashregexwarnings .= "Sequence: $identifier, $seqlength bp\n";
    print "\nSequence: $identifier, $seqlength bp\n";

    # replace all ?s in the sequence with Ns
    $hashofsequences{$identifier} =~ s/\?/N/g;

    # verify that the sequence ends with a complete codon - the sequence length should be evenly divisible by 3
    if ($seqlength%3 != 0) {
      # write error to error log file
      $warnings .= "Warning: Sequence ($identifier, $seqlength bp) is not in the first reading frame or does not end with a complete codon\n";
      print "Warning: Sequence ($identifier, $seqlength bp) is not in the first reading frame or does not end with a complete codon\n";
    }

    # Jump through the sequence and check each codon
    for($position = 0; $position < $seqlength; $position += 3) {
      # Extract the codon at the current position as upper case characters from the sequence
      $codon1 = uc substr($hashofsequences{$identifier}, $position, 3);
      # Compare the codon to our lists (in a hash) and find the degenerate codon with which to replace it
      if (exists $degencodons1by1{$codon1}) {
        $degensequences{$identifier} .= "$degencodons1by1{$codon1}";
        $hashregexwarnings .= "hash: $codon1 replaced with $degencodons1by1{$codon1} at position " . ($position + 1) . "\n";
      } else {
        if ($codon1 =~ m/([^-][RYMKSWHBVDN]|[RYMKSWHBVDN][^-])([^-])/i) {
          # if nt1 and/or nt2 are/is degenerate, the codon becomes NNN
          $hashregexwarnings .= "Warning: $codon1 replaced with ";
          print "Warning: $codon1 replaced with ";
          $codon1 =~ s/([^-][RYMKSWHBVDN]|[RYMKSWHBVDN][^-])(.)/NNN/i;
          $degensequences{$identifier} .= "$codon1";
          $hashregexwarnings .= "$codon1 at position " . ($position + 1) . "\n";
          print "$codon1 at position " . ($position + 1) . "\n";
        } else {
            $warnings .= "Warning: Codon ($codon1, " . ($position + 1) . ") not changed\n";
            print "Warning: Codon ($codon1, " . ($position + 1) . ") not changed\n";
            $degensequences{$identifier} .= "$codon1";
        }
      }
    }
  }
} # END degenerate


sub create_output () { # write results and errors to files
  # New FASTA file with degenerate sequences
  open (DEGEN, ">Degen_$filename.fasta") or die "An error occured while attempting to open (W) Degen_$filename/Degen_$filename.fasta: $!";
  foreach (sort keys %degensequences) {
    print DEGEN ">$_\n$degensequences{$_}\n";
  }
  close(DEGEN);

  # New NEXUS file with degenerate sequences
  open (NEXUS, ">Degen_$filename.nex") or die "An error occured while attempting to open (W) Degen_$filename/Degen_$filename.nex: $!";
  print NEXUS "#NEXUS\n";
  print NEXUS "Begin Data;\n";
  print NEXUS "  Dimensions ntax=", scalar (keys %degensequences), " nchar=$seqlength;\n";
  print NEXUS "  Format datatype=dna gap=-;\n";
  print NEXUS "  Matrix\n";
  foreach (sort keys %degensequences) {
    print NEXUS "    $_\t$degensequences{$_}\n";
  }
  print NEXUS "  ;\n";
  print NEXUS "END;\n";
  close(NEXUS);

  # New warning output file
  open(WARN, ">Warnings_$filename.txt") or die "An error occured while attempting to open (W) Degen_$filename/Warnings_$filename.txt: $!";
  print WARN $warnings;
  close (WARN);

  # New hash/regex output file
  open(HASHRE, ">HashRegEx_$filename.txt") or die "An error occured while attempting to open (W) Degen_$filename/HashRegEx_$filename.txt: $!";
  print HASHRE $hashregexwarnings;
  close (HASHRE);
} # END create_output


sub help_text () {
  my $error = $_[0];
  print "\n DEGENERACY CODING $version \t(c) Andreas Zwick & April Hussey 2008-2012\n";
  print "\************************\n";
  if ( ! $error eq "") {
    print "\n$error\n";
  }
  print "\nCorrect usage:\n\n";
  print "\tDegen.pl <input file> [--table=#]\n\n";
  print "REQUIRED:\n";
  print "  'input file': Sequence file in interleaved or non-interleaved FASTA, FLAT or NEXUS format.\n";
  print "\nOPTIONAL:\n";
  print "  -t= or --table=#: NCBI codon table number (only 1-6, 9-14 are supported) or 'S', 'Z' or 'SZ'\n";
  print "                   (Standard genetic code [1] is default)\n";
  exit;
}#END help_text
