                             Degen v1.4
 
           Andreas Zwick, April Hussey & Jerome C. Regier
      State Museum of Natural History Stuttgart, Germany (AZ)
               University of Maryland, USA (AH, JCR)
 
     Comments or questions about this script should be sent to
              Andreas Zwick at degen1@phylotools.com.
 
README file for Degen.pl, version 1.4
Last updated: 5 JAN 2013
Latest version available at: http://www.phylotools.com
 
 
*********************************************************************
*  Please acknowledge the use of this script in your publications   *
*  by citing:                                                       *
*                                                                   *
*  Zwick, A., Regier, J.C. & Zwickl, D.J. (2012).                   *
*     "Resolving Discrepancy between Nucleotides and Amino Acids    *
*     in Deep-Level Arthropod Phylogenomics: Differentiating        *
*     Serine Codons in 21-Amino-Acid Models".                       *
*     PLoS ONE 7(11): e47450.                                       *
*                                                                   *
*  Regier, J.C., Shultz, J.W., Zwick, A., Hussey, A., Ball, B.,     *
*     Wetzer, R. Martin, J.W. & Cunningham, C.W. (2010).            *
*     "Arthropod relationships revealed by phylogenomic analysis    *
*     of nuclear protein-coding sequences". Nature 463: 1079-1083.  *
*                                                                   *
*********************************************************************
 
 
CONTENTS:
1) Introduction
2) Requirements and installation
3) Usage
4) How it works
5) Versions
6) History of development
7) License
8) Acknowledgments
 
 
1) INTRODUCTION
================
The reconstruction of evolutionary relationships between ancient
lineages using highly divergent, protein-coding DNA sequence data can
be hampered by nucleotide compositional heterogeneity. Current
versions of software frequently used to infer molecular phylogenies
(e.g., MrBayes, RAxML, Garli, PAUP*) do not account for compositional
heterogeneity across taxa, with divergences from homogeneity
resulting in signals that can conflict, but occasionally concur, with
the phylogenetic signal inferred from the sequence of nucleotides.
     The PERL script "Degen" aids in the phylogenetic analysis of
highly divergent DNA sequence data by greatly reducing nucleotide
compositional heterogeneity between taxa. The key observation that
"Degen" exploits is that most compositional heterogeneity resides in
sites that undergo synonymous change. "Degen" operates by
degenerating nucleotides at all sites that can potentially undergo
synonymous change in any and all pairwise comparisons of sequences in
the data matrix, thereby making synonymous change largely invisible
and reducing compositional heterogeneity but leaving the inference of
nonsynonymous change largely intact. The "Degen" script fully
degenerates all codons that encode single amino acids by substituting
one of the four standard nucleotides (A,C,G,T) with IUPAC ambiguity
codes (M,R,W,S,Y,K,V,H,D,B,N) that allow for all possible synonymous
change for that amino acid. In the current version of "Degen"
(v1.4), the "standard", various mitochondrial, plastid and other 
genetic codes for nuclear, protein-coding genes in animals, plants 
and other organisms is implemented. A web-based script that allows 
users to directly transform their input data is available at 
http://www.phylotools.com.
     For more background information and details, see section 4, 
"HOW IT WORKS".
 
 
2) REQUIREMENTS & INSTALLATION
===============================
Data requirements:
      - nucleotide sequences of protein-coding genes that conform 
        to a supported genetic code [codes 1-6 and 9-14]:
        http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
      - sequences have to be in first open reading frame and
        consist only of complete codons
        [sequence begins with 1st and ends with 3rd codon position]
      - any indels, represented by dashes in the data matrix, must
          be triplets or their multiples and in-frame, i.e., located
          between nt3 and nt1 of different codons
      - sequence data are in FASTA, FLAT or NEXUS file format
      - typically, but not necessarily, the data file is a
        multi-sequence alignment
 
System requirements:
      - any operating system (e.g., Linux/UNIX, Mac OS-X, Windows)
        with a functional installation of a PERL language interpreter
        (e.g., http://www.activestate.com/activeperl/);
        type "perl -V" in a shell to check for a PERL installation
 
There is no need to install the script other than to copy it to a
directory of your choice. To use the script directly from any
directory, place the script in a location that is included in your
PATH variable or adjust the PATH variable accordingly.
 
 
3) USAGE
=========
The script expects only a single command-line parameter, namely the
data input file:
       Degen_v1_4.pl <datafilename>
 
Optionally, a genetic code (1-6 and 9-14 of NCBI) or alternative 
encoding for the "standard genetic code" (S, Z or SZ) can be 
specified (default is 1, the "standard genetic code"):
       Degen_v1_4.pl <datafilename> [-t= OR --table=##]
 
It is best called with your PERL interpreter, but the actual command
might vary depending on your operating system.
 
An example of how to use the script correctly under Linux with the 
invertebrate mitochondrial code:
       perl ./Degen_v1_4.pl mydata.nex --table=5
 
The output doesn't overwrite the original input file, but consists 
of four files in a single folder ["Degen_datafilename"]. These files
are the new, degenerated matrix in FASTA and NEXUS format
["Degen_datafilename.fasta" and "Degen_datafilename.nex"], a list 
of each and every nucleotide transformation 
["HashRegEx_datafilename.txt"] and a list of sequence lengths and 
potentially problematic nucleotide triplets such as termination 
codons, out-of-frame indels and unexpected characters 
["Warnings_datafilename.txt"].
 
 
4) HOW IT WORKS
================
The "Degen" script reads individual DNA sequences as strings of 
codons, in which there are three sequential nucleotides per codon
(nt1 nt2 nt3). It then replaces every codon with a fully
"degenerated" codon, using IUPAC nomenclature of polymorphic
nucleotides (e.g., C+T = Y) for those nucleotides that can be
variable, yet encode for the same amino acid. Hence, whenever there
are multiple codons that encode the same amino acid, the original
nucleotides of the input sequence are expanded to match all such
codons.
     At nt2, nucleotides are left unaltered, since no synonymous
differences based on single nucleotide substitutions are possible at
that position. Most of the expansions occur at the highly variable
nt3, for example, GGG --> GGN (glycine), ATT --> ATH (isoleucine),
GAT --> GAY (aspartic acid) but ATG (methionine) remains as is.
In addition to degenerating synonymous differences, the script also
modifies already polymorphic codons that encode more than one amino
acid (typically not commonplace). Triplets that encode leucine +
phenylalanine are converted to YTN. Triplets that encode arginine +
serine2 are converted to MGN. Triplets that encode histidine +
glutamine, asparagine + lysine and aspartic acid + glutamic acid are
converted to CAN, AAN and GAN, respectively. All other nucleotide
triplets that encode such non-synonymous polymorphisms are converted
to NNN.
     Any indels in the sequence, represented by (multiples of) three
dashes (---), are not modified by the script. In contrast, missing
data, represented by question marks, are converted to N's.
 
Standard genetic code [1]:
For leucine and arginine codons, of which there are six each, nt3 is
converted to "N". For nt1, however, there is no fully satisfactory
manner to degenerate leucine- and arginine-encoding nucleotides using
single, informative tokens without also being slightly misinformative
or without eliminating synonymous and nonsynonymous differences. We
have explored three approaches (see below) and find that the first
one, which is the implemented standard, is the most effective.
     This first approach ("Degen1") is to convert all leucine-encoding 
codons to YTN, and all arginine-encoding sequences to MGN. By this 
approach the conversion of, e.g., TTG and CTG (both leucine) to YTN 
encodes not only leucine, but also phenylalanine (TTT and TTC). 
Likewise, the conversion of, e.g., AGA and CGA (both arginine) to MGN 
encodes not only arginine, but also serine2 (AGT and AGC = AGY, as 
compared to serine1, which is encoded by TCN).
     A second approach, which would eliminate the just-mentioned
"phenylalanine" and "serine2" miscoding problems, is to leave the nt1
codings of leucine and arginine unaltered. For example, by converting
leucine TTG --> TTR (leucine) and CTG --> CTN (leucine), but this
would retain synonymous differences at nt1 across all leucine codons
(T <--> C). The equivalent would be the case for arginine (AGG -->
AGR and CGG --> CGN), which retains synonymous differences at nt1
across all arginine codons (A <--> C).
     A third approach is to convert all nt1 characters that encode
either for leucine / arginine or for phenylalanine / serine2 to "N",
but this results in a substantial loss of information.
 
The effect of "Degen" coding is to minimize, but not necessarily to
eliminate, the detection of synonymous change through the
introduction of an analytical criterion for data set modification.
This criterion is that the input data matrix be degenerated such that
when inferring change in a parsimony framework between any and all
pairs of terminal sequences, synonymous change must become invisible
(undetectable) while nonsynonymous change be left largely intact.
There are two sequelae of particular note. Firstly, because
"Degen1" degenerates leucine to YTN and arginine to MGN, it actually
reinterprets leucine as part leucine and part phenylalanine, and
arginine as part arginine and part serine2 (see above). However, this
affects only selected nt1 characters, as it is nucleotides and not
codons that are analyzed. 
We note that currently available software cannot analyze a "Degen" 
matrix under a codon model due to the use of ambiguity codes. Secondly, 
because "Degen" degenerates selected sites, it changes the estimate of 
the average base composition in analysis software, such as PAUP* and Garli.
 
 
5) VERSIONS
============
The latest version of the "Degen1" script and this README file can be
downloaded at http://www.phylotools.com
 
Version 1.4:  5 JAN 2013
Version 1.3: 19 NOV 2009
Version 1.2: 14 MAY 2009
Version 1.1: 10 SEP 2008
Version 1.0: 26 AUG 2008
 
 
6) HISTORY OF DEVELOPMENT
==========================
In 2008, Jerome Regier, Andreas Zwick and April Hussey discussed how 
more nonsynonymous change than that yielded by "LeuArg" analysis 
might be extracted from a data matrix without increasing synonymous 
change. The resulting concept was "degen1". April Hussey wrote the 
first script and then modified it in 2009 (Degen1.pl, version 1.2). 
It was first published in: Regier, J.C., Shultz, J.W., Zwick, A., 
Hussey, A. Ball, B., Wetzer, R., Martin, J.W. Cunningham, C. (2010).
Arthropod relationships revealed by phylogenomic anlaysis of nuclear 
protein-coding sequences. Nature. In Press.
 
Andreas Zwick updated and expanded the script (version 1.3 - 1.4), adding 
support for alternative genetic codes, alternative encodings of the 
"standard genetic codes", Nexus files and interleaved format. He also 
makes the script and its web service available at:
http://www.phylotools.com
 
 
7) LICENSE
===========
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version. This program is distributed in the
hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.
 
 
8) ACKNOWLEDGMENTS
====================
Development of the "degen1" script (through February, 2010) was 
funded by a grant from the U.S. National Science Foundation 
(Biocomplexity in the Environment: Genome-Enabled Environmental 
Science and Engineering program, Award no. DEB-0120635; Assembling 
the Tree of Life program, Award no. 0531626).
