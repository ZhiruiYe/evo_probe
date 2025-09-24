      seqfile = final_align_fixed.phy  * pruned PHYLIP alignment (amino acids) *
     treefile = Rename_newick_tree.nwk  * pruned Newick tree with branch lengths *
      outfile = mlc 

        noisy = 3
      verbose = 2
      runmode = 0                        * user tree *

        seqtype = 2                      * 2 = amino acids *
       aaRatefile = /yezhirui/evo_probe/paml-4.10.9-linux-x86_64/dat/jones.dat            * JTT matrix *
           clock = 0
          model = 2                      * empirical AA model *
         NSsites = 0

     fix_blength = 1                     * use given branch lengths *
      cleandata = 0                      * keep gaps/missing as '?' *

     fix_alpha = 0                       * estimate gamma shape (alpha) *
          alpha = 1.0                    * initial value *
          ncatG = 4                      * JTT+Gamma(4) *

       getSE = 0
 RateAncestor = 1                        * reconstruct ancestors; output to 'rst' *
