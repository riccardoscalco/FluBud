Flubud calculates the Rohlin distance between alphabetical strings, according
to a certain procedure of partitioning.

For details, see:
Burioni R , Scalco R , Casartelli M (2011) Rohlin Distance and the Evolution of
Influenza A Virus: Weak Attractors and Precursors. PLoS ONE 6(12): e27924.
doi:10.1371/journal.pone.0027924


Flags
-----

--s : The accepted difference between labels (optional, default to 1).
      
      Here there are some examples with the string "aababbac":
      s=1 --> [[1,2],[3],[4],[5,6],[7],[8]]
      s=2 --> [[1,2,4],[3,5,6],[7],[8]]
      s=3 --> [[1,2,4,7],[3,5,6],[8]]

--fasta : The fasta file to read as input (required).
          The file has to have the following format:
          - sequence must be of the same length;
          - the description line of the sequence must start with '>' and cannot
            start a new line (one line description);
          - empty lines are not allowed (note that some times there is an empty
            line at the end of the file);

          Here there is an example:

------------8<--------------------------------------------------------
>AEO08998 A/Ankara/01/2010 2010/12/ HA
MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLG
KCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFINYEELREQLSSVSSFERFEIFPKTSSW
PNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQ
NADAYVFVGTSKYSKKFKPEIAVRPKVRDQEGRMDYYWTLVEPGDKITFEATGNLLVPRYAFAMERNAGS
GIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAI
AGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDKITNKVNSVIEKMNTQFTAVGKEFNHLEKR
IENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDN
TCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLIVSLGAISFWMCSNGSL
QCRICI
>AEM92338 A/District of Columbia/WRAIR0306/2010 2010/12/28 HA
MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLG
KCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSW
PNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLNKSYINDKGKEVLVLWGIHHPSTTADQQSLYQ
NADAYVFVGTSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGS
GIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAI
AGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDKITNKVNSVIEKMNTQFTAVGKEFNHLEKR
IENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRNQLKNNAKEIGNGCFEFYHKCDN
TCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSL
QCRICI
>AEF28931 A/Guangdong/002/2011 2011/01/19 HA
MKAILVVMLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLG
KCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSW
PNHDSNKGVTAACPHAGAKGFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTTADQQSLYQ
NADTYVFVGTSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGS
GFIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAI
AGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDKITNKVNSVIEKMNTQFTAVGKEFNHLEKR
IENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRNQLKNNAKEIGNGCFEFYHKCDN
TCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSL
QCRICI
------------8<--------------------------------------------------------


Output
------

The output is the Rohlin matrix. The matrix is symmetric with zeros
on the diagonal, therefore only the up-diagonal elements are printed.

This is the output of the above example:

[5.545177444479549,8.841014310483871]
[3.295836866004322]

Where the first element of the first list corresponds to the element M[1,2]
in the matrix, the second element of the first list corresponds to M[1,3] and
the first element of the second list corresponds to M[2,3].

Usage
-----

The following command calculates the Rohlin matrix for the set of sequences
listed on the "FASTA.fa" file, with s=10 as maximum difference between labels.

$ ./flubud --s=10 --fasta="FASTA.fa" > out &

Observations
------------

1) There is no percentage indicator for the actual job to be done,
so for who of you that are impatient it is suggested to use the standard unix
commands "grep" and "wc": first, find the total number of sequences in your file
FASTA.fa with

$ grep ">" FASTA.fa | wc

(the first number is the number of sequences), then do

$ wc out

and you will see the first number growing to the number of sequences.


2) The measure used in the calculation of entropy is not normalized to
one. It is instead normalized to the length of the sequences. However, the
normalization to one is restored simply by dividing the matrix (the output)
by the length of the sequences.

To do
-----

- let the user to define the measure to be used in the calculation of entropy.

- introduce the reduction process for arbitrary partitions.
